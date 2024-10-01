!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_estimate_condition

  private
  public :: hecmw_estimate_condition_CG
  public :: hecmw_estimate_condition_GMRES

contains

  subroutine hecmw_estimate_condition_CG(ITER, D, E)
    use hecmw_util
    implicit none
    integer(kind=kint), intent(in) :: ITER
    real(kind=kreal), intent(in) :: D(:), E(:)

#ifdef HECMW_WITH_LAPACK
    character(len=1) :: JOBZ, range
    ! character(len=1) :: COMPZ
    real(kind=kreal) :: VL, VU, ABSTOL, Z(1,1)
    integer(kind=kint) :: N, IL, IU, M, LDZ=1, ISUPPZ(1)
    integer(kind=kint) :: LWORK, LIWORK, INFO
    real(kind=kreal), allocatable :: W(:), WORK(:)
    integer(kind=kint), allocatable :: IWORK(:)
    real(kind=kreal), allocatable :: D1(:), E1(:)
    integer(kind=kint) :: i

    if (ITER <= 1) return

    ! copy D, E
    allocate(D1(ITER),E1(ITER))
    do i=1,ITER-1
      D1(i) = D(i)
      E1(i) = E(i)
    enddo
    D1(ITER) = D(ITER)


    !!
    !! dstegr version (faster than dsteqr)
    !!

    ! prepare arguments for calling dstegr
    JOBZ='N'
    range='A'
    N=ITER
    allocate(W(ITER))
    ! estimate optimal LWORK and LIWORK
    LWORK=-1
    LIWORK=-1
    allocate(WORK(1),IWORK(1))
    call dstegr(JOBZ,range,N,D1,E1,VL,VU,IL,IU,ABSTOL, &
      M,W,Z,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: dstegr returned with INFO=',INFO
      return
    endif
    ! calculate eigenvalues
    LWORK=WORK(1)
    LIWORK=IWORK(1)
    deallocate(WORK,IWORK)
    allocate(WORK(LWORK),IWORK(LIWORK))
    call dstegr(JOBZ,range,N,D1,E1,VL,VU,IL,IU,ABSTOL, &
      M,W,Z,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: dstegr returned with INFO=',INFO
      return
    endif
    write(*,'("emin=",1pe13.6,", emax=",1pe13.6,", emax/emin=",1pe13.6)') &
      W(1),W(N),W(N)/W(1)
    deallocate(WORK,IWORK)
    deallocate(W)


    ! !!
    ! !! dsteqr version
    ! !!

    ! ! prepare arguments for calling dsteqr
    ! COMPZ='N'
    ! N=ITER
    ! allocate(WORK(1))
    ! ! calculate eigenvalues
    ! call dsteqr(COMPZ,N,D1,E1,Z,LDZ,WORK,INFO)
    ! if (INFO /= 0) then
    !   write(*,*) 'ERROR: dsteqr returned with INFO=',INFO
    !   return
    ! endif
    ! write(*,'("emin=",1pe13.6,", emax=",1pe13.6,", emax/emin=",1pe13.6)') &
      !      D1(1),D1(N),D1(N)/D1(1)
    ! deallocate(WORK)


    deallocate(D1,E1)
#endif
  end subroutine hecmw_estimate_condition_CG


  subroutine hecmw_estimate_condition_GMRES(I, H)
    use hecmw_util
    implicit none
    integer(kind=kint), intent(in) :: I
    real(kind=kreal), intent(in) :: H(:,:)

#ifdef HECMW_WITH_LAPACK
    ! character(len=1) :: JOBU, JOBVT
    character(len=1) :: JOBZ
    integer(kind=kint) :: N, LDH, LDZ=1, LWORK, INFO
    real(kind=kreal), allocatable :: WR(:), WORK(:), H1(:,:)
    integer(kind=kint), allocatable :: IWORK(:)
    real(kind=kreal) :: Z(1,1)
    integer(kind=kint) :: j, k

    if (I == 0) return

    ! copy H
    N=I
    allocate(H1(N+1,N))
    do j = 1, N
      do k = 1, j+1
        H1(k,j) = H(k,j)
      enddo
      do k = j+2, N
        H1(k,j) = 0.d0
      enddo
    enddo
    LDH=N+1
    allocate(WR(N))


    ! !!
    ! !! dgesvd version
    ! !!

    ! ! arguments for calling dgesvd
    ! JOBU='N'
    ! JOBVT='N'
    ! ! estimate optimal LWORK
    ! allocate(WORK(1))
    ! LWORK=-1
    ! call dgesvd(JOBU,JOBVT,N,N,H1,LDH,WR,Z,LDZ,Z,LDZ,WORK,LWORK,INFO)
    ! if (INFO /= 0) then
    !   write(*,*) 'ERROR: dgesvd returned with INFO=',INFO
    !   return
    ! endif
    ! ! calculate singular values
    ! LWORK=WORK(1)
    ! deallocate(WORK)
    ! allocate(WORK(LWORK))
    ! call dgesvd(JOBU,JOBVT,N,N,H1,LDH,WR,Z,LDZ,Z,LDZ,WORK,LWORK,INFO)
    ! if (INFO /= 0) then
    !   write(*,*) 'ERROR: dgesvd returned with INFO=',INFO
    !   return
    ! endif
    ! deallocate(WORK)


    !!
    !! dgesdd version (faster but need more workspace than dgesvd)
    !!

    ! arguments for calling dgesdd
    JOBZ='N'
    allocate(IWORK(8*N))
    ! estimate optimal LWORK
    allocate(WORK(1))
    LWORK=-1
    call dgesdd(JOBZ,N,N,H1,LDH,WR,Z,LDZ,Z,LDZ,WORK,LWORK,IWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: dgesdd returned with INFO=',INFO
      return
    endif
    ! calculate singular values
    LWORK=WORK(1)
    deallocate(WORK)
    allocate(WORK(LWORK))
    call dgesdd(JOBZ,N,N,H1,LDH,WR,Z,LDZ,Z,LDZ,WORK,LWORK,IWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: dgesdd returned with INFO=',INFO
      return
    endif
    deallocate(WORK)
    deallocate(IWORK)


    write(*,'("emin=",1pe13.6,", emax=",1pe13.6,", emax/emin=",1pe13.6)') &
      WR(N), WR(1), WR(1)/WR(N)

    deallocate(WR)
    deallocate(H1)
#endif
  end subroutine hecmw_estimate_condition_GMRES

end module hecmw_estimate_condition
