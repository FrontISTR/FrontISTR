!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_SOLVE_DIRECT is a program for the matrix direct solver
!         SOLVER=DIRECT
!----------------------------------------------------------------------
module HECMW_SOLVER_DIRECT
  use HECMW_UTIL
  implicit none

  private
  public :: HECMW_SOLVE_DIRECT

  real(kind=kreal), parameter :: RMIn = 4.941D-300

  integer(kind=kint), parameter :: IDBg_ini = 0
  integer(kind=kint), parameter :: IDBg_sym = 0
  integer(kind=kint), parameter :: IDBg_num = 0

  type cholesky_factor
    integer(kind=kint) :: LEN_colno
    integer(kind=kint) :: NSTop
    integer(kind=kint) :: STAge
    integer(kind=kint) :: NEQns
    integer(kind=kint) :: NDEg
    integer(kind=kint) :: NTTbr
    integer(kind=kint) :: ISYm
    integer(kind=kint) :: LEN_dsln
    integer(kind=kint), pointer :: JCOl(:)
    integer(kind=kint), pointer :: IROw(:)
    integer(kind=kint), pointer :: IPErm(:)
    integer(kind=kint), pointer :: INVp(:)
    integer(kind=kint), pointer :: PARent(:)
    integer(kind=kint), pointer :: NCH(:)
    integer(kind=kint), pointer :: XLNzr(:)
    integer(kind=kint), pointer :: COLno(:)
    real(kind=kreal), pointer :: DIAg(:)
    real(kind=kreal), pointer :: ZLN(:)
    real(kind=kreal), pointer :: DSLn(:)
    !*Allocation variables
    integer(kind=kint) :: ialoc
    integer(kind=kint) :: raloc
  end type cholesky_factor

contains
  !----------------------------------------------------------------------
  !> @brief HECMW_SOLVE_DIRECT is a program for the matrix solver
  !     and is also available for performance tests
  !     to solve Ax=b for x, nozero pattern and values must be give.
  !   arrays
  !     jcol     column entries
  !     irow     row entries
  !     val      its value
  !     v        interface array used through matini, staijx, nufctx,
  !              nusolx
  !----------------------------------------------------------------------
  subroutine HECMW_SOLVE_DIRECT(hecMESH,hecMAT,Ifmsg)
    use HECMW_MATRIX_ASS
    use HECMW_MATRIX_DUMP
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    type (HECMWST_MATRIX), intent(inout)::hecMAT
    integer(kind=kint), intent(in):: Ifmsg
    !------
    type (cholesky_factor), save :: FCT
    !------
    integer(kind=kint):: i98
    integer(kind=kint):: i97
    integer(kind=kint):: timelog
    integer(kind=kint):: iterlog
    integer(kind=kint):: ordering
    integer(kind=kint):: loglevel
    integer(kind=kint):: ir
    integer(kind=kint):: i
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5

    ir = 0

    timelog = hecMAT%IARRAY(22)
    iterlog = hecMAT%IARRAY(21)
    ordering = hecMAT%IARRAY(41)
    loglevel = max(timelog,iterlog)

    call HECMW_MAT_DUMP(hecMAT,hecMESH)

    call PTIME(t1)
    t2 = t1

    !*EHM HECMW June 7 2004
    i98 = hecMAT%IARRAY(98)
    if ( hecMAT%IARRAY(98)==1 ) then
      !* Interface to symbolic factorization
      call SETIJ(hecMESH,hecMAT,FCT)

      !* Symbolic factorization
      call MATINI(FCT,ordering,loglevel,ir)
      hecMAT%IARRAY(98) = 0

      if ( loglevel > 0  ) write (*,*) "[DIRECT]: symbolic fct done"
    endif
    call PTIME(t2)
    t3 = t2

    i97 = hecMAT%IARRAY(97)
    if ( hecMAT%IARRAY(97)==1 ) then
      !* Interface to numeric factorization
      call NUFORM(hecMESH,hecMAT,FCT,ir)
      call PTIME(t3)

      !* Numeric factorization
      call NUFCT0(FCT,ir)
      hecMAT%IARRAY(97) = 0

      if ( loglevel > 0 ) write (*,*) "[DIRECT]: numeric fct done"

      !*Memory Details
      if ( loglevel > 1 ) then
        write (*,*) '*-----------------------------------*'
        write (*,*) '|   Direct  Solver  Memory  Usage   |'
        write (*,*) '*-----------------------------------*'
        write (*,*) 'INTEGER memory: ', real(FCT%IALoc*4)/real(1048576), 'MB'
        write (*,*) 'REAL*8  memory: ', real(FCT%RALoc*8)/real(1048576), 'MB'
        write (*,*) 'TOTAL   memory: ', real((FCT%RALoc*2+FCT%IALoc)*4)/real(1048576), 'MB'
        write (*,*) '*-----------------------------------*'
      endif
    endif
    call PTIME(t4)

    !* Finalize
    !*  Errors 1
    if ( i98/=0 .and. i98/=1 ) then
      write (Ifmsg,*) 'ERROR in symb. fact. flag: Should be 1 or 0'
      stop 'ERROR in symb. fact. flag: Should be 1 or 0'
    endif
    if ( i97/=0 .and. i97/=1 ) then
      write (Ifmsg,*) 'ERROR in numer. fact. flag: Should be 1 or 0'
      stop 'ERROR in numer. fact. flag: Should be 1 or 0'
    endif
    if ( i98==1 .and. i97==0 ) then
      write (Ifmsg,*) 'WARNING: Numeric factorization not performed!'
      stop 'WARNING: Numeric factorization not performed! Solve will not be performed'
    endif
    !*  Errors 2
    if ( ir/=0 ) then
      write (Ifmsg,*) 'ERROR in nufct0. ir = ', ir
      stop
    endif

    !* Solve
    do i=1,hecMAT%NP*hecMESH%n_dof
      hecMAT%X(i) = hecMAT%B(i)
    end do
    call NUSOL0(hecMAT%X,FCT,ir)
    call PTIME(t5)
    !* Errors 4
    if ( ir/=0 ) then
      write (Ifmsg,*) 'error in nusol0. irr = ', ir
      stop
    endif

    if ( timelog > 1 ) then
      write(*,*) 'sym fct time : ',t2-t1
      write(*,*) 'nuform time  : ',t3-t2
      write(*,*) 'num fct time : ',t4-t3
      write(*,*) 'solve time   : ',t5-t4
    elseif ( timelog > 0 ) then
      write(*,*) 'setup time : ',t4-t1
      write(*,*) 'solve time : ',t5-t4
    endif

    call HECMW_MAT_DUMP_SOLUTION(hecMAT)
  end subroutine HECMW_SOLVE_DIRECT

  !======================================================================!
  !> @brief PTIME
  !======================================================================!
  subroutine PTIME(Cputim)
    implicit none
    !------
    real(kind=kreal), intent(out):: Cputim
    !------
    ! cpu time by hour
    Cputim = HECMW_WTIME()
  end subroutine PTIME

  !======================================================================!
  !> @brief SETIJ
  !  sets NEQns, NDEg, NTTbr, ISYm, JCOl, IROw
  !======================================================================!
  subroutine SETIJ(hecMESH,hecMAT,FCT)
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    type (HECMWST_MATRIX), intent(in)::hecMAT
    type (cholesky_factor), intent(inout) :: FCT
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ierr
    integer(kind=kint):: j
    integer(kind=kint):: k
    integer(kind=kint):: kk
    integer(kind=kint):: ntotal
    integer(kind=kint):: numnp
    integer(kind=kint):: ndof
    integer(kind=kint):: ndof2

    numnp = hecMAT%NP
    ndof = hecMESH%N_DOF
    ntotal = numnp*ndof

    !*NUFACT variables
    FCT%NEQns = numnp
    FCT%NDEg = ndof
    FCT%NTTbr = hecMAT%NP + hecMAT%NPL
    !+hecMAT%NPU if unsymmetric
    FCT%ISYm = 0

    !*Allocations
    allocate (FCT%IROw(FCT%NTTbr),stat=ierr)
    if ( ierr/=0 ) stop "Allocation error: irow"
    allocate (FCT%JCOl(FCT%NTTbr),stat=ierr)
    if ( ierr/=0 ) stop "Allocation error: jcol"

    kk = 0
    ndof2 = ndof*ndof
    do j = 1, numnp
      !*Diagonal
      kk = kk + 1
      FCT%IROw(kk) = j
      FCT%JCOl(kk) = j
      !*Lower
      do k = hecMAT%INDEXL(j-1) + 1, hecMAT%INDEXL(j)
        i = hecMAT%ITEML(k)
        kk = kk + 1
        FCT%IROw(kk) = j
        FCT%JCOl(kk) = i
      enddo
    enddo
  end subroutine SETIJ

  !======================================================================!
  !> @brief MATINI initializes storage for sparse matrix solver.
  !     this routine is used for both symmetric and asymmetric matrices
  !     and must be called once at the beginning
  !    (i)
  !        neqns     number of unknowns
  !        nttbr     number of non0s, pattern of non-zero elements are
  !                  given like following.
  !                  nonz(A)={(i,j);i=irow(l),j=jcol(l); 1<= l <= nttbr}
  !        irow
  !        jcol      to define non-zero pattern
  !        lenv      length of the array v (iv)
  !    (o)
  !        iv        communication array. v is the original name
  !        ir        return code
  !                              =0    normal
  !                              =-1   non positive index
  !                              =1    too big index
  !                              =10   insufficient storage
  !        contents of iv
  !               pointers 1 &zpiv(1)  2 &iperm(1)  3 &invp(1)
  !                        4 &parent(1)5 &nch(1)    6 &xlnzr(1)
  !                        7 &colno(1) 8 &diag(1)   9 &zln(1)
  !                       10 &dsln(1)
  !
  !               scalars 21 len(colno)  22 nstop     23 stage
  !                       24 neqns       25 len(iv) 26 len(dsln)
  !                       27 total
  !        stage   10  after initialization
  !                20  building up matrix
  !                30  after LU decomposition
  !                40  after solving
  !
  !  sets LEN_colno, NSTop, LEN_dsln, IPErm, INVp, PARent, NCH, XLNzr, COLno
  !       IALoc, RALoc, STAge
  !======================================================================!
  subroutine MATINI(FCT,ordering,loglevel,Ir)
    use hecmw_ordering
    implicit none
    !------
    type (cholesky_factor), intent(inout) :: FCT
    integer(kind=kint), intent(in):: ordering
    integer(kind=kint), intent(in):: loglevel
    integer(kind=kint), intent(out):: Ir
    !------
    !*Work arrays
    integer(kind=kint), allocatable :: zpiv(:)
    integer(kind=kint), allocatable :: jcpt(:)
    integer(kind=kint), allocatable :: jcolno(:)
    integer(kind=kint), allocatable :: ia(:)
    integer(kind=kint), allocatable :: ja(:)
    integer(kind=kint), allocatable :: quarent(:)
    integer(kind=kint), allocatable :: btree(:)
    integer(kind=kint), allocatable :: xleaf(:)
    integer(kind=kint), allocatable :: leaf(:)
    integer(kind=kint):: ir1
    integer(kind=kint):: irr
    integer(kind=kint):: izz
    integer(kind=kint):: izz0
    integer(kind=kint):: lncol
    integer(kind=kint):: neqnsz
    integer(kind=kint):: ierror

    izz0 = 0

    Ir = 0

    !*Initialize allocation measure variables
    FCT%IALoc = 0
    FCT%RALoc = 0
    !
    !  set z pivot
    !
    allocate (ZPIv(FCT%NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, zpiv: SUB. matini"
    call ZPIVOT(FCT%NEQns,neqnsz,FCT%NTTbr,FCT%JCOl,FCT%IROw,ZPIv,ir1)
    if ( ir1/=0 ) then
      Ir = ir1
      return
    endif

    !
    !  build jcpt,jcolno
    !
    allocate (JCPt(2*FCT%NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, jcpt: SUB. matini"
    allocate (JCOlno(2*FCT%NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, jcolno: SUB. matini"
    call STSMAT(FCT%NEQns,FCT%NTTbr,FCT%IROw,FCT%JCOl,JCPt,JCOlno)
    !
    !  build ia,ja
    !
    allocate (IA(FCT%NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, ia: SUB. matini"
    allocate (JA(2*FCT%NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, ja: SUB. matini"
    call STIAJA(FCT%NEQns,IA,JA,JCPt,JCOlno)

    !*Deallocation of work array
    deallocate (JCPt)
    deallocate (JCOlno)
    !
    !  get permutation vector iperm,invp
    !
    allocate (FCT%IPErm(FCT%NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, iperm: SUB. matini"
    allocate (FCT%INVp(FCT%NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, invp: SUB. matini"
    allocate (Quarent(FCT%NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, quarent: SUB. matini"
    call hecmw_ordering_GEN(neqnsz,FCT%NTTbr,IA,JA,FCT%IPErm,FCT%INVp,ordering,loglevel)
    do
      !   build up the parent vector
      !   parent vector will be saved in Quarent for a while
      call GENPAQ(IA,JA,FCT%INVp,FCT%IPErm,Quarent,FCT%NEQns)
      !
      !   build up the binary tree
      !
      allocate (BTRee(2*(FCT%NEQns+1)),stat=IERror)
      if ( IERror/=0 ) stop "ALLOCATION ERROR, btree: SUB. matini"
      call GENBTQ(FCT%INVp,Quarent,BTRee,ZPIv,izz,FCT%NEQns)
      !
      !   rotate the binary tree to avoid a zero pivot
      !
      if ( izz==0 ) then
        !
        !   post ordering
        !
        allocate (FCT%PARent(FCT%NEQns),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, parent: SUB. matini.f"
        allocate (FCT%NCH(FCT%NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, nch: SUB. matini.f"
        call POSORD(FCT%PARent,BTRee,FCT%INVp,FCT%IPErm,FCT%NCH,FCT%NEQns,Quarent)
        !
        !   generate skeleton graph
        !
        allocate (XLEaf(FCT%NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, xleaf: SUB. matini.f"
        allocate (LEAf(FCT%NTTbr),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, leaf: SUB. matini.f"
        call GNLEAF(IA,JA,FCT%INVp,FCT%IPErm,FCT%NCH,XLEaf,LEAf,FCT%NEQns)
        call FORPAR(FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop)
        !*Deallocation of work arrays
        deallocate (IA)
        deallocate (JA)
        deallocate (Quarent)
        deallocate (ZPIv)
        !
        !   build up xlnzr,colno  (this is the symbolic fct.)
        !
        allocate (FCT%XLNzr(FCT%NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, xlnzr: SUB. matini.f"
        call PRE_GNCLNO(FCT%PARent,XLEaf,LEAf,FCT%XLNzr,FCT%NEQns,FCT%NSTop,lncol,ir1)
        allocate (FCT%COLno(lncol),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, colno: SUB. matini.f"
        call GNCLNO(FCT%PARent,XLEaf,LEAf,FCT%XLNzr,FCT%COLno,FCT%NEQns,FCT%NSTop,lncol,ir1)
        !*Deallocate work arrays
        deallocate (XLEaf)
        deallocate (LEAf)
        deallocate (BTRee)
        FCT%LEN_dsln = (FCT%NEQns-FCT%NSTop+1)*(FCT%NEQns-FCT%NSTop)/2

        !Scalar assignments
        FCT%LEN_colno = lncol

        FCT%STAge = 10
        FCT%IALoc = 5*FCT%NEQns + lncol + 1
        exit
      else
        if ( izz0==0 ) izz0 = izz
        if ( izz0/=izz ) then
          call BRINGU(ZPIv,FCT%IPErm,FCT%INVp,Quarent,izz,FCT%NEQns,IRR)
        else
          call ROTATE(IA,JA,FCT%INVp,FCT%IPErm,Quarent,BTRee,izz,FCT%NEQns,IRR)
        endif
      endif
    enddo
  end subroutine MATINI

  !======================================================================!
  !> @brief ZPIVOT
  !======================================================================!
  subroutine ZPIVOT(Neqns,Neqnsz,Nttbr,Jcol,Irow,Zpiv,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nttbr
    integer(kind=kint), intent(in):: Jcol(:)
    integer(kind=kint), intent(in):: Irow(:)
    integer(kind=kint), intent(out):: Ir
    integer(kind=kint), intent(out):: Neqnsz
    integer(kind=kint), intent(out):: Zpiv(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: l

    Ir = 0
    Zpiv(1:Neqns) = 1

    loop1: do
      do l = 1, Nttbr
        i = Irow(l)
        j = Jcol(l)
        if ( i<=0 .or. j<=0 ) then
          Ir = -1
          exit loop1
        elseif ( i>Neqns .or. j>Neqns ) then
          Ir = 1
          exit loop1
        endif
        if ( i==j ) Zpiv(i) = 0
      enddo
      do i = Neqns, 1, -1
        if ( Zpiv(i)==0 ) then
          Neqnsz = i
          exit
        endif
      enddo
      exit
    enddo loop1
    if ( IDBg_ini/=0 ) write (6,"(20I3)") (Zpiv(i),i=1,Neqns)
  end subroutine ZPIVOT

  !======================================================================!
  !> @brief STSMAT
  !======================================================================!
  subroutine STSMAT(Neqns,Nttbr,Irow,Jcol,Jcpt,Jcolno)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nttbr
    integer(kind=kint), intent(in):: Irow(:)
    integer(kind=kint), intent(in):: Jcol(:)
    integer(kind=kint), intent(out):: Jcpt(:)
    integer(kind=kint), intent(out):: Jcolno(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: l
    integer(kind=kint):: locr
    logical:: found

    Jcpt(1:2*Nttbr) = 0
    Jcolno(1:2*Nttbr) = 0
    do i = 1, Neqns
      Jcpt(i) = i + Neqns
      Jcolno(i+Neqns) = i
    enddo

    k = 2*Neqns

    loop1: do l = 1, Nttbr
      i = Irow(l)
      j = Jcol(l)
      if ( i/=j ) then
        joc = Jcpt(i)
        locr = i
        found = .false.
        do while ( joc/=0 )
          if ( Jcolno(joc)==j ) cycle loop1
          if ( Jcolno(joc)>j ) then
            k = k + 1
            Jcpt(locr) = k
            Jcpt(k) = joc
            Jcolno(k) = j
            found = .true.
            exit
          endif
          locr = joc
          joc = Jcpt(joc)
        enddo
        if (.not. found) then
          k = k + 1
          Jcpt(locr) = k
          Jcolno(k) = j
        endif

        joc = Jcpt(j)
        locr = j
        do while ( joc/=0 )
          if ( Jcolno(joc)==i ) cycle loop1
          if ( Jcolno(joc)>i ) then
            k = k + 1
            Jcpt(locr) = k
            Jcpt(k) = joc
            Jcolno(k) = i
            cycle loop1
          endif
          locr = joc
          joc = Jcpt(joc)
        enddo
        k = k + 1
        Jcpt(locr) = k
        Jcolno(k) = i
      endif
    enddo loop1
    if ( IDBg_ini/=0 ) then
      write (6,*) 'jcolno'
      write (6,"(10I7)") (Jcolno(i),i=1,k)
      write (6,*) 'jcpt'
      write (6,"(10I7)") (Jcpt(i),i=1,k)
    endif
  end subroutine STSMAT

  !======================================================================!
  !> @brief STIAJA routine sets an non-zero entry  of the matrix.
  !      (asymmetric version)
  !======================================================================!
  subroutine STIAJA(Neqns,Ia,Ja,Jcpt,Jcolno)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Jcpt(:)
    integer(kind=kint), intent(in):: Jcolno(:)
    integer(kind=kint), intent(out):: Ia(:)
    integer(kind=kint), intent(out):: Ja(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ii
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: l

    Ia(1) = 1
    l = 0
    do k = 1, Neqns
      joc = Jcpt(k)
      do while ( joc/=0 )
        ii = Jcolno(joc)
        if ( ii/=k ) then
          l = l + 1
          Ja(l) = ii
        endif
        joc = Jcpt(joc)
      enddo
      Ia(k+1) = l + 1
    enddo
    if ( IDBg_ini/=0 ) then
      write (6,*) 'ia '
      write (6,"(10I7)") (Ia(i),i=1,Neqns)
      write (6,*) 'ja '
      write (6,"(10I7)") (Ja(i),i=1,Ia(Neqns+1))
    endif
  end subroutine STIAJA

  !======================================================================!
  !> @brief GENPAQ
  !======================================================================!
  subroutine GENPAQ(Xadj,Adjncy,Invp,Iperm,Parent,Neqns)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Iperm(:)
    integer(kind=kint), intent(out):: Parent(:)
    !------
    integer(kind=kint), allocatable:: Ancstr(:)
    integer(kind=kint):: i
    integer(kind=kint):: ip
    integer(kind=kint):: it
    integer(kind=kint):: k
    integer(kind=kint):: l
    integer(kind=kint):: ierror

    allocate (Ancstr(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, ancstr: SUB. genpaq"

    do i = 1, Neqns
      Parent(i) = 0
      Ancstr(i) = 0
      ip = Iperm(i)
      loop1: do k = Xadj(ip), Xadj(ip+1) - 1
        l = Invp(Adjncy(k))
        if ( l<i ) then
          do while ( Ancstr(l)/=0 )
            if ( Ancstr(l)==i ) cycle loop1
            it = Ancstr(l)
            Ancstr(l) = i
            l = it
          enddo
          Ancstr(l) = i
          Parent(l) = i
        endif
      enddo loop1
    enddo
    do i = 1, Neqns
      if ( Parent(i)==0 ) Parent(i) = Neqns + 1
    enddo
    Parent(Neqns+1) = 0

    if ( IDBg_sym/=0 ) then
      write (6,"(' parent')")
      write (6,"(2I6)") (i,Parent(i),i=1,Neqns)
    endif

    deallocate(Ancstr)
  end subroutine GENPAQ

  !======================================================================!
  !> @brief GENBTQ
  !======================================================================!
  subroutine GENBTQ(Invp,Parent,Btree,Zpiv,Izz,Neqns)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Zpiv(:)
    integer(kind=kint), intent(out):: Btree(2,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ib
    integer(kind=kint):: inext
    integer(kind=kint):: ip
    integer(kind=kint):: Izz

    Btree(1,1:Neqns + 1) = 0
    Btree(2,1:Neqns + 1) = 0
    do i = 1, Neqns + 1
      ip = Parent(i)
      if ( ip>0 ) then
        ib = Btree(1,ip)
        if ( ib==0 ) then
          Btree(1,ip) = i
        else
          do
            inext = Btree(2,ib)
            if ( inext==0 ) then
              Btree(2,ib) = i
            else
              ib = inext
              cycle
            endif
            exit
          enddo
        endif
      endif
    enddo
    !
    ! find zeropivot
    !
    Izz = 0
    do i = 1, Neqns
      if ( Zpiv(i)/=0 ) then
        if ( Btree(1,Invp(i))==0 ) then
          Izz = i
          exit
        endif
      endif
    enddo

    if ( IDBg_sym/=0 ) then
      write (6,"(' binary tree')")
      write (6,"(i6,'(',2I6,')')") (i,Btree(1,i),Btree(2,i),i=1,Neqns)
      write (6,"(' the first zero pivot is ',i4)") Izz
    endif
  end subroutine GENBTQ

  !======================================================================!
  !> @brief POSORD
  !======================================================================!
  subroutine POSORD(Parent,Btree,Invp,Iperm,Nch,Neqns,Qarent)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Btree(2,*)
    integer(kind=kint), intent(in):: Qarent(:)
    integer(kind=kint), intent(out):: Parent(:)
    integer(kind=kint), intent(inout):: Invp(:)
    integer(kind=kint), intent(out):: Iperm(:)
    integer(kind=kint), intent(out):: Nch(:)
    !------
    integer(kind=kint), allocatable:: Pordr(:)
    integer(kind=kint), allocatable:: Iw(:)
    integer(kind=kint), allocatable:: Mch(:)
    integer(kind=kint):: i
    integer(kind=kint):: ii
    integer(kind=kint):: invpos
    integer(kind=kint):: ipinv
    integer(kind=kint):: joc
    integer(kind=kint):: l
    integer(kind=kint):: locc
    integer(kind=kint):: locp
    integer(kind=kint):: ierror

    allocate (Pordr(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, pordr: SUB. posord"
    allocate (Iw(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, iw: SUB. posord"
    allocate (Mch(0:Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, mch: SUB. posord"

    Mch(1:Neqns) = 0
    Pordr(1:Neqns) = 0
    l = 1
    locc = Neqns + 1
    do
      joc = locc
      locc = Btree(1,joc)
      if ( locc==0 ) then
        locp = Qarent(joc)
        Mch(locp) = Mch(locp) + 1
        do
          Pordr(joc) = l
          if ( l>=Neqns ) then
            do i = 1, Neqns
              ipinv = Pordr(Invp(i))
              Invp(i) = ipinv
              Iperm(ipinv) = i
              Iw(Pordr(i)) = i
            enddo
            do i = 1, Neqns
              invpos = Iw(i)
              Nch(i) = Mch(invpos)
              ii = Qarent(invpos)
              if ( ii>0 .and. ii<=Neqns ) then
                Parent(i) = Pordr(ii)
              else
                Parent(i) = Qarent(invpos)
              endif
            enddo
            if ( IDBg_sym/=0 ) then
              write (6,"(' post order')")
              write (6,"(10I6)") (Pordr(i),i=1,Neqns)
              write (6,"(/' invp ')")
              write (6,"(/' parent')")
              write (6,"(10I6)") (Parent(i),i=1,Neqns)
              write (6,"(10I6)") (Invp(i),i=1,Neqns)
              write (6,"(/' iperm ')")
              write (6,"(10I6)") (Iperm(i),i=1,Neqns)
              write (6,"(' nch')")
              write (6,"(10I6)") (Nch(i),i=1,Neqns)
            endif
            return
          else
            l = l + 1
            locc = Btree(2,joc)
            if ( locc/=0 ) exit
            joc = Qarent(joc)
            locp = Qarent(joc)
            Mch(locp) = Mch(locp) + Mch(joc) + 1
          endif
        enddo
      endif
    enddo

    deallocate (Pordr)
    deallocate (Iw)
    deallocate (Mch)
  end subroutine POSORD

  !======================================================================!
  !> @brief GNLEAF
  !======================================================================!
  subroutine GNLEAF(Xadj,Adjncy,Invp,Iperm,Nch,Xleaf,Leaf,Neqns)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: Nch(:)
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Iperm(:)
    integer(kind=kint), intent(out):: Xleaf(:)
    integer(kind=kint), intent(out):: Leaf(:)
    !------
    integer(kind=kint), allocatable:: Adjncp(:)
    integer(kind=kint):: Lnleaf
    integer(kind=kint):: i
    integer(kind=kint):: ik
    integer(kind=kint):: ip
    integer(kind=kint):: iq
    integer(kind=kint):: istart
    integer(kind=kint):: k
    integer(kind=kint):: l
    integer(kind=kint):: lc
    integer(kind=kint):: lc1
    integer(kind=kint):: m
    integer(kind=kint):: ierror

    allocate (Adjncp(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, adjncp: SUB. gnleaf"

    l = 1
    ik = 0
    istart = 0
    do i = 1, Neqns
      Xleaf(i) = l
      ip = Iperm(i)
      do k = Xadj(ip), Xadj(ip+1) - 1
        iq = Invp(Adjncy(k))
        if ( iq<i ) then
          ik = ik + 1
          Adjncp(ik) = iq
        endif
      enddo
      m = ik - istart
      if ( m/=0 ) then
        call QQSORT(Adjncp(istart+1:),m)
        lc1 = Adjncp(istart+1)
        if ( lc1<i ) then
          Leaf(l) = lc1
          l = l + 1
          do k = istart + 2, ik
            lc = Adjncp(k)
            if ( lc1<lc-Nch(lc) ) then
              Leaf(l) = lc
              l = l + 1
            endif

            lc1 = lc
          enddo
          ik = 1
          istart = ik
        endif
      endif
    enddo
    Xleaf(Neqns+1) = l
    Lnleaf = l - 1

    if ( IDBg_sym/=0 ) then
      write (6,"(' xleaf')")
      write (6,"(10I6)") (Xleaf(i),i=1,Neqns+1)
      write (6,"(' leaf (len = ',i6,')')") Lnleaf
      write (6,"(10I6)") (Leaf(i),i=1,Lnleaf)
    endif

    deallocate (Adjncp)
    return
  end subroutine GNLEAF

  !======================================================================!
  !> @brief QQSORT
  ! sort in increasing order up to i
  !     iw   array
  !     ik   number of input/output
  !     i    deal with numbers less than this numberi
  !======================================================================!
  subroutine QQSORT(Iw,Ik)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ik
    integer(kind=kint), intent(inout):: Iw(:)
    !------
    integer(kind=kint):: itemp
    integer(kind=kint):: l
    integer(kind=kint):: m

    if ( Ik<=1 ) return
    do l = 1, Ik - 1
      do m = l + 1, Ik
        if ( Iw(l)>=Iw(m) ) then
          itemp = Iw(l)
          Iw(l) = Iw(m)
          Iw(m) = itemp
        endif
      enddo
    enddo
  end subroutine QQSORT

  !======================================================================!
  !> @brief FORPAR
  !======================================================================!
  subroutine FORPAR(Neqns,Parent,Nch,Nstop)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Nstop
    integer(kind=kint), intent(out):: Nch(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: idens
    integer(kind=kint):: ii

    Nch(1:Neqns) = 0
    Nch(Neqns+1) = 0
    do i = 1, Neqns
      ii = Parent(i)
      Nch(ii) = Nch(ii) + 1
    enddo
    do i = Neqns, 1, -1
      if ( Nch(i)/=1 ) exit
    enddo

    idens = 0
    if ( idens==1 ) then
      Nstop = i
    else
      Nstop = Neqns + 1
    endif
  end subroutine FORPAR

  !======================================================================!
  !> @brief PRE_GNCLNO
  !======================================================================!
  subroutine PRE_GNCLNO(Parent,Xleaf,Leaf,Xlnzr,Neqns,Nstop,Lncol,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(in):: Xleaf(:)
    integer(kind=kint), intent(in):: Leaf(:)
    integer(kind=kint), intent(out):: Lncol
    integer(kind=kint), intent(out):: Xlnzr(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: Ir
    integer(kind=kint):: j
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    integer(kind=kint):: nc
    integer(kind=kint):: nxleaf

    nc = 0
    Ir = 0
    l = 1
    loop1: do i = 1, Neqns
      Xlnzr(i) = l
      ks = Xleaf(i)
      ke = Xleaf(i+1) - 1
      if ( ke>=ks ) then
        nxleaf = Leaf(ks)
        do k = ks, ke - 1
          j = nxleaf
          nxleaf = Leaf(k+1)
          do while ( j<nxleaf )
            if ( j>=Nstop ) cycle loop1
            l = l + 1
            j = Parent(j)
          enddo
        enddo
        j = Leaf(ke)
        do while ( j<Nstop )
          if ( j>=i .or. j==0 ) exit
          l = l + 1
          j = Parent(j)
        enddo
      endif
    enddo loop1
    Xlnzr(Neqns+1) = l
    Lncol = l - 1
  end subroutine PRE_GNCLNO

  !======================================================================!
  !> @brief GNCLNO
  !======================================================================!
  subroutine GNCLNO(Parent,Xleaf,Leaf,Xlnzr,Colno,Neqns,Nstop,Lncol,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(in):: Xleaf(:)
    integer(kind=kint), intent(in):: Leaf(:)
    integer(kind=kint), intent(out):: Ir
    integer(kind=kint), intent(out):: Lncol
    integer(kind=kint), intent(out):: Xlnzr(:)
    integer(kind=kint), intent(out):: Colno(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    integer(kind=kint):: nc
    integer(kind=kint):: nxleaf

    nc = 0
    Ir = 0
    l = 1
    loop1: do i = 1, Neqns
      Xlnzr(i) = l
      ks = Xleaf(i)
      ke = Xleaf(i+1) - 1
      if ( ke>=ks ) then
        nxleaf = Leaf(ks)
        do k = ks, ke - 1
          j = nxleaf
          nxleaf = Leaf(k+1)
          do while ( j<nxleaf )
            if ( j>=Nstop ) cycle loop1
            Colno(l) = j
            l = l + 1
            j = Parent(j)
          enddo
        enddo
        j = Leaf(ke)
        do while ( j<Nstop )
          if ( j>=i .or. j==0 ) exit
          Colno(l) = j
          l = l + 1
          j = Parent(j)
        enddo
      endif
    enddo loop1
    Xlnzr(Neqns+1) = l
    Lncol = l - 1

    if ( IDBg_sym/=0 ) then
      write (6,"(' xlnzr')")
      write (6,"(' colno (lncol =',i10,')')") Lncol
      do k = 1, Neqns
        write (6,"(/' row = ',i6)") k
        write (6,"(10I4)") (Colno(i),i=Xlnzr(k),Xlnzr(k+1)-1)
      enddo
    endif
  end subroutine GNCLNO

  !======================================================================!
  !> @brief BRINGU  brings up zero pivots from bottom of the elimination tree to higher nodes
  !      irr = 0     complete
  !          = 1     impossible
  !======================================================================!
  subroutine BRINGU(Zpiv,Iperm,Invp,Parent,Izz,Neqns,Irr)
    implicit none
    !------
    integer(kind=kint), intent(in):: Izz
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Zpiv(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Irr
    integer(kind=kint), intent(inout):: Iperm(:)
    integer(kind=kint), intent(inout):: Invp(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ib
    integer(kind=kint):: ib0
    integer(kind=kint):: ibp
    integer(kind=kint):: idbg
    integer(kind=kint):: izzp

    idbg = 0
    Irr = 0
    ib0 = Invp(Izz)
    ib = ib0
    do while ( ib>0 )
      ibp = Parent(ib)
      izzp = Iperm(ibp)
      if ( Zpiv(izzp)==0 ) then
        Invp(Izz) = ibp
        Invp(izzp) = ib0
        Iperm(ibp) = Izz
        Iperm(ib0) = izzp
        if ( idbg/=0 ) then
          do i = 1, Neqns
            if ( Invp(Iperm(i))/=i .or. Iperm(Invp(i))/=i) then
              write (6,*) 'permutation error'
              stop
            endif
          enddo
          return
        endif
        return
      else
        ib = ibp
      endif
    enddo
    Irr = 1
  end subroutine BRINGU

  !======================================================================!
  !> @brief ROTATE
  ! irr return code irr=0 node izz is not a bottom node
  !                     irr=1          is a bottom node then rotation is
  !                                    performed
  !======================================================================!
  subroutine ROTATE(Xadj,Adjncy,Invp,Iperm,Parent,Btree,Izz,Neqns,Irr)
    implicit none
    !------
    integer(kind=kint), intent(in):: Izz
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(in):: Btree(2,*)
    integer(kind=kint), intent(out):: Irr
    integer(kind=kint), intent(inout):: Invp(:)
    integer(kind=kint), intent(inout):: Iperm(:)
    !------
    integer(kind=kint), allocatable:: Anc(:)
    integer(kind=kint), allocatable:: Adjt(:)
    integer(kind=kint):: i
    integer(kind=kint):: iy
    integer(kind=kint):: izzz
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: l
    integer(kind=kint):: ll
    integer(kind=kint):: locc
    integer(kind=kint):: nanc
    integer(kind=kint):: ierror

    allocate (Anc(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, anc: SUB. rotate"
    allocate (Adjt(Neqns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, adjt: SUB. rotate"

    if ( Izz==0 ) then
      Irr = 0
      return
    endif
    izzz = Invp(Izz)

    if ( Btree(1,izzz)/=0 ) Irr = 0
    Irr = 1
    !
    !  ancestors of izzz
    !
    nanc = 0
    joc = izzz
    do
      nanc = nanc + 1
      Anc(nanc) = joc
      joc = Parent(joc)
      if ( joc==0 ) then
        !
        !  to find the eligible node from ancestors of izz
        !
        l = 1
        exit
      endif
    enddo

    loop1: do
      Adjt(1:Neqns) = 0
      locc = Anc(l)
      do
        joc = locc
        locc = Btree(1,joc)
        if ( locc==0 ) then
          do
            Adjt(Invp(Adjncy(Xadj(Iperm(joc)):Xadj(Iperm(joc)+1) - 1))) = 1
            if ( joc>=Anc(l) ) then
              do ll = l + 1, nanc
                if ( Adjt(Anc(ll))==0 ) then
                  l = l + 1
                  cycle loop1
                endif
              enddo
              if ( l==1 ) then
                !
                ! izz can be numbered last
                !
                k = 0
                do i = 1, Neqns
                  if ( i/=izzz ) then
                    k = k + 1
                    Invp(Iperm(i)) = k
                  endif
                enddo
                Invp(Iperm(izzz)) = Neqns
              else
                !
                !  anc(l-1) is the eligible node
                !
                ! (1) number the node not in Ancestor(iy)
                iy = Anc(l-1)
                Adjt(1:Neqns) = 0
                Adjt(Anc(l:nanc)) = 1
                k = 0
                do ll = 1, Neqns
                  if ( Adjt(ll)==0 ) then
                    k = k + 1
                    Invp(Iperm(ll)) = k
                  endif
                enddo
                ! (2) followed by nodes in Ancestor(iy)-Adj(T(iy))
                Adjt(1:Neqns) = 0
                locc = iy
                loop2: do
                  joc = locc
                  locc = Btree(1,joc)
                  if ( locc==0 ) then
                    do
                      Adjt(Invp(Adjncy(Xadj(Iperm(joc)):Xadj(Iperm(joc)+1)-1))) = 1
                      if ( joc>=iy ) then
                        do ll = l, nanc
                          if ( Adjt(Anc(ll))==0 ) then
                            k = k + 1
                            Invp(Iperm(Anc(ll))) = k
                          endif
                        enddo
                        ! (3) and finally number the node in Adj(t(iy))
                        do ll = l, nanc
                          if ( Adjt(Anc(ll))/=0 ) then
                            k = k + 1
                            Invp(Iperm(Anc(ll))) = k
                          endif
                        enddo
                        exit loop2
                      else
                        locc = Btree(2,joc)
                        if ( locc/=0 ) exit
                        joc = Parent(joc)
                      endif
                    enddo
                  endif
                enddo loop2
              endif
              !
              ! set iperm
              !
              do i = 1, Neqns
                Iperm(Invp(i)) = i
              enddo

              if ( IDBg_sym/=0 ) write (6,"(10I6)") (Invp(i),i=1,Neqns)
              return
            else
              locc = Btree(2,joc)
              if ( locc/=0 ) exit
              joc = Parent(joc)
            endif
          enddo
        endif
      enddo
      exit
    enddo loop1

    deallocate (Anc)
    deallocate (Adjt)
  end subroutine ROTATE

  !======================================================================!
  !> @brief NUFORM
  !  sets DIAg, ZLN, DSLn
  !  updates RALoc, STAge
  !======================================================================!
  subroutine NUFORM(hecMESH,hecMAT,FCT,Ir)
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    type (HECMWST_MATRIX), intent(in)::hecMAT
    type (cholesky_factor), intent(inout):: FCT
    integer(kind=kint), intent(out):: Ir
    !------
    integer(kind=kint):: i
    integer(kind=kint):: idbg
    integer(kind=kint):: ierr
    integer(kind=kint):: j
    integer(kind=kint):: k
    integer(kind=kint):: kk
    integer(kind=kint):: ntotal
    integer(kind=kint):: numnp
    integer(kind=kint):: ndof
    integer(kind=kint):: ndof2
    real(kind=kreal), allocatable :: val(:)

    idbg = 0
    numnp = hecMAT%NP
    ndof = hecMESH%N_DOF
    ntotal = numnp*ndof

    !*Allocations
    allocate (val(FCT%NDEg*FCT%NDEg),stat=ierr)
    if ( ierr/=0 ) stop "Allocation error:val"
    if ( IDBg_num/= 0 ) write (6,*) "nuform:stage = ", FCT%STAge
    kk = 0
    ndof2 = ndof*ndof

    do j = 1, numnp
      !*Diagonal
      kk = kk + 1
      call VLCPY(val,hecMAT%D(ndof2*(j-1)+1:ndof2*j),ndof)
      call STAIJ1(0,j,j,val,FCT,Ir)

      do i = 1, ndof
        if ( val((i-1)*ndof+i)<=0 ) write (idbg,*) 'j,j,val:', j, i, val((i-1)*ndof+i)
      enddo

      !*Lower
      do k = hecMAT%INDEXL(j-1) + 1, hecMAT%INDEXL(j)
        i = hecMAT%ITEML(k)
        kk = kk + 1
        call VLCPY(val,hecMAT%AL(ndof2*(k-1)+1:ndof2*k),ndof)
        call STAIJ1(0,j,i,val,FCT,Ir)
      enddo
    enddo

    deallocate (val)
  end subroutine NUFORM

  !======================================================================!
  !> @brief VLCPY
  !======================================================================!
  subroutine VLCPY(A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: B(N,N)
    real(kind=kreal), intent(out):: A(N,N)
    !------
    integer(kind=kint):: i

    do i = 1, N
      A(1:N,i) = B(i,1:N)
    enddo
  end subroutine VLCPY

  !======================================================================!
  !> @brief STAIJ1 routine sets an non-zero entry  of the matrix.
  !      (symmetric version)
  !      (i)
  !          isw      =0    set the value
  !                   =1    add the value
  !          i        row entry
  !          j        column entry
  !          aij      value
  !      (o)
  !          iv       communication array
  !  sets DIAg, ZLN, DSLn
  !  updates RALoc, STAge
  !======================================================================!
  subroutine STAIJ1(Isw,I,J,Aij,FCT,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: I
    integer(kind=kint), intent(in):: Isw
    integer(kind=kint), intent(in):: J
    real(kind=kreal), intent(in):: Aij(:)
    type (cholesky_factor), intent(inout):: FCT
    integer(kind=kint), intent(out):: Ir
    !------
    integer(kind=kint):: ndeg2
    integer(kind=kint):: ndeg2l
    integer(kind=kint):: ierror

    Ir = 0
    ndeg2 = FCT%NDEg*FCT%NDEg
    ndeg2l = FCT%NDEg*(FCT%NDEg+1)/2
    if ( FCT%STAge==30 ) write (6,*) 'warning a matrix was build up '//'but never solved.'
    if ( FCT%STAge==10 ) then
      allocate (FCT%DIAg(FCT%NEQns*ndeg2l),stat=IERror)
      if ( IERror/=0 ) stop "Allocation error diag"
      FCT%RALoc = FCT%RALoc + FCT%NEQns*ndeg2l

      allocate (FCT%ZLN(FCT%LEN_colno*ndeg2),stat=IERror)
      if ( IERror/=0 ) stop "Allocation error zln"
      FCT%RALoc = FCT%RALoc + FCT%LEN_colno*ndeg2

      allocate (FCT%DSLn(FCT%LEN_dsln*ndeg2),stat=IERror)
      if ( IERror/=0 ) stop "Allocation error dsln"
      FCT%RALoc = FCT%RALoc + FCT%LEN_dsln*ndeg2
    endif
    if ( FCT%STAge/=20 ) then
      !
      ! for diagonal
      !
      FCT%DIAg = 0.
      !
      ! for lower triangle
      !
      FCT%ZLN = 0.
      !
      ! for dense window
      !
      FCT%DSLn = 0.

      FCT%STAge = 20
    endif
    !         Print *,'********Set Stage 20 *********'
    !
    if ( FCT%NDEg<=2 ) then
      call ADDR0(Isw,I,J,Aij,FCT%INVp,FCT%XLNzr,FCT%COLno,FCT%DIAg,FCT%ZLN,FCT%DSLn,FCT%NSTop,ndeg2,ndeg2l,Ir)
    elseif ( FCT%NDEg==3 ) then
      call ADDR3(I,J,Aij,FCT%INVp,FCT%XLNzr,FCT%COLno,FCT%DIAg,FCT%ZLN,FCT%DSLn,FCT%NSTop,Ir)
    else
      call ADDRX(I,J,Aij,FCT%INVp,FCT%XLNzr,FCT%COLno,FCT%DIAg,FCT%ZLN,FCT%DSLn,FCT%NSTop,FCT%NDEg,ndeg2l,Ir)
    endif
  end subroutine STAIJ1

  !======================================================================!
  !> @brief ADDR0
  !======================================================================!
  subroutine ADDR0(Isw,I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ndeg2,Ndeg2l,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: I
    integer(kind=kint), intent(in):: Isw
    integer(kind=kint), intent(in):: J
    integer(kind=kint), intent(in):: Ndeg2
    integer(kind=kint), intent(in):: Ndeg2l
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Aij(Ndeg2)
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Zln(Ndeg2,*)
    real(kind=kreal), intent(inout):: Diag(Ndeg2l,*)
    real(kind=kreal), intent(inout):: Dsln(Ndeg2,*)
    !------
    integer(kind=kint):: i0
    integer(kind=kint):: ii
    integer(kind=kint):: itrans
    integer(kind=kint):: j0
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint), parameter:: idbg = 0

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    if ( idbg/=0 ) write (6,*) ii, jj, Aij
    if ( ii==jj ) then
      if ( Ndeg2==1 ) then
        if ( Isw==0 ) then
          Diag(1,ii) = Aij(1)
        else
          Diag(1,ii) = Diag(1,ii) + Aij(1)
        endif
      elseif ( Ndeg2==4 ) then
        if ( Isw==0 ) then
          Diag(1,ii) = Aij(1)
          Diag(2,ii) = Aij(2)
          Diag(3,ii) = Aij(4)
        else
          Diag(1,ii) = Diag(1,ii) + Aij(1)
          Diag(2,ii) = Diag(2,ii) + Aij(2)
          Diag(3,ii) = Diag(3,ii) + Aij(4)
        endif
      endif
      return
    endif
    itrans = 0
    if ( jj>ii ) then
      k = jj
      jj = ii
      ii = k
      itrans = 1
    endif
    if ( jj>=Nstop ) then
      i0 = ii-Nstop
      j0 = jj-Nstop + 1
      k = i0*(i0-1)/2 + j0
      if ( Ndeg2==1 ) then
        Dsln(1,k) = Aij(1)
        return
      elseif ( Ndeg2==4 ) then
        if ( itrans==0 ) then
          Dsln(1:Ndeg2,k) = Aij(1:Ndeg2)
        else
          Dsln(1,k) = Aij(1)
          Dsln(2,k) = Aij(3)
          Dsln(3,k) = Aij(2)
          Dsln(4,k) = Aij(4)
        endif
        return
      endif
    endif
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    do k = ks, ke
      if ( Colno(k)==jj ) then
        if ( Isw==0 ) then
          if ( Ndeg2==1 ) then
            Zln(1,k) = Aij(1)
          elseif ( Ndeg2==4 ) then
            if ( itrans==0 ) then
              Zln(1:Ndeg2,k) = Aij(1:Ndeg2)
            else
              Zln(1,k) = Aij(1)
              Zln(2,k) = Aij(3)
              Zln(3,k) = Aij(2)
              Zln(4,k) = Aij(4)
            endif
          endif
        elseif ( Ndeg2==1 ) then
          Zln(1,k) = Zln(1,k) + Aij(1)
        elseif ( Ndeg2==4 ) then
          if ( itrans==0 ) then
            Zln(1:Ndeg2,k) = Zln(1:Ndeg2,k) + Aij(1:Ndeg2)
          else
            Zln(1,k) = Zln(1,k) + Aij(1)
            Zln(2,k) = Zln(2,k) + Aij(3)
            Zln(3,k) = Zln(3,k) + Aij(2)
            Zln(4,k) = Zln(4,k) + Aij(4)
          endif
        endif
        return
      endif
    enddo
    Ir = 20
  end subroutine ADDR0

  !======================================================================!
  !> @brief ADDR3
  !======================================================================!
  subroutine ADDR3(I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: I
    integer(kind=kint), intent(in):: J
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Aij(9)
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Zln(9,*)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(inout):: Dsln(9,*)
    !------
    integer(kind=kint):: i0
    integer(kind=kint):: ii
    integer(kind=kint):: itrans
    integer(kind=kint):: j0
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    integer(kind=kint), parameter:: idbg = 0
    integer(kind=kint), parameter:: ndeg2 = 9
    integer(kind=kint), parameter:: ndeg2l = 6

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    if ( idbg/=0 ) write (6,*) ii, jj, Aij
    if ( ii==jj ) then
      Diag(1,ii) = Aij(1)
      Diag(2,ii) = Aij(2)
      Diag(3,ii) = Aij(5)
      Diag(4,ii) = Aij(3)
      Diag(5,ii) = Aij(6)
      Diag(6,ii) = Aij(9)
      return
    endif
    itrans = 0
    if ( jj>ii ) then
      k = jj
      jj = ii
      ii = k
      itrans = 1
    endif
    if ( jj>=Nstop ) then
      i0 = ii - Nstop
      j0 = jj - Nstop + 1
      k = i0*(i0-1)/2 + j0
      if ( itrans==0 ) then
        Dsln(1:ndeg2,k) = Aij(1:ndeg2)
      else
        Dsln(1,k) = Aij(1)
        Dsln(2,k) = Aij(4)
        Dsln(3,k) = Aij(7)
        Dsln(4,k) = Aij(2)
        Dsln(5,k) = Aij(5)
        Dsln(6,k) = Aij(8)
        Dsln(7,k) = Aij(3)
        Dsln(8,k) = Aij(6)
        Dsln(9,k) = Aij(9)
      endif
      return
    endif
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    do k = ks, ke
      if ( Colno(k)==jj ) then
        if ( itrans==0 ) then
          do l = 1, ndeg2
            Zln(l,k) = Aij(l)
          enddo
        else
          Zln(1,k) = Aij(1)
          Zln(2,k) = Aij(4)
          Zln(3,k) = Aij(7)
          Zln(4,k) = Aij(2)
          Zln(5,k) = Aij(5)
          Zln(6,k) = Aij(8)
          Zln(7,k) = Aij(3)
          Zln(8,k) = Aij(6)
          Zln(9,k) = Aij(9)
        endif
        return
      endif
    enddo
    Ir = 20
  end subroutine ADDR3

  !======================================================================!
  !> @brief ADDRX
  !======================================================================!
  subroutine ADDRX(I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ndeg,Ndeg2l,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: I
    integer(kind=kint), intent(in):: J
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndeg2l
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Invp(:)
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Aij(Ndeg,Ndeg)
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Zln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Diag(Ndeg2l,*)
    real(kind=kreal), intent(inout):: Dsln(Ndeg,Ndeg,*)
    !------
    integer(kind=kint):: i0
    integer(kind=kint):: ii
    integer(kind=kint):: itrans
    integer(kind=kint):: j0
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    integer(kind=kint):: m
    integer(kind=kint):: n
    integer(kind=kint), parameter:: idbg = 0

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    if ( idbg/=0 ) write (6,*) ii, jj, Aij
    if ( ii==jj ) then
      l = 0
      do n = 1, Ndeg
        do m = 1, n
          l = l + 1
          Diag(l,ii) = Aij(n,m)
        enddo
      enddo
      return
    endif
    itrans = 0
    if ( jj>ii ) then
      k = jj
      jj = ii
      ii = k
      itrans = 1
    endif
    if ( jj>=Nstop ) then
      i0 = ii - Nstop
      j0 = jj - Nstop + 1
      k = i0*(i0-1)/2 + j0
      if ( itrans==0 ) then
        do m = 1, Ndeg
          Dsln(1:Ndeg,m,k) = Aij(1:Ndeg,m)
        enddo
      else
        do m = 1, Ndeg
          Dsln(1:Ndeg,m,k) = Aij(m,1:Ndeg)
        enddo
      endif
      return
    endif
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    do k = ks, ke
      if ( Colno(k)==jj ) then
        if ( itrans==0 ) then
          do m = 1, Ndeg
            Zln(1:Ndeg,m,k) = Aij(1:Ndeg,m)
          enddo
        else
          do m = 1, Ndeg
            Zln(1:Ndeg,m,k) = Aij(m,1:Ndeg)
          enddo
        endif
        return
      endif
    enddo
    Ir = 20
  end subroutine ADDRX

  !======================================================================!
  !> @brief NUFCT0 performs Cholesky factorization
  !          if(iv(22).eq.0)    normal type
  !          if(iv(22).gt.0)    code generation type
  !  updates NCH, DIAg, ZLN, DSLn, STAge
  !======================================================================!
  subroutine NUFCT0(FCT,Ir)
    implicit none
    !------
    type (cholesky_factor), intent(inout):: FCT
    integer(kind=kint), intent(out):: Ir
    !------
    integer(kind=kint):: irr
    integer(kind=kint):: ndeg2
    integer(kind=kint):: ndegl
    integer(kind=kint), allocatable :: INDx(:)
    real(kind=kreal), allocatable :: TEMp(:)

    if ( FCT%STAge/=20 ) then
      print *, '*********Setting Stage 40!*********'
      Ir = 40
      return
    else
      Ir = 0
    endif

    allocate (TEMp(FCT%NDEg*FCT%NDEg*FCT%NEQns),stat=IRR)
    if ( IRR/=0 ) then
      write (*,*) '##Error : Not enough memory'
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
      !stop
    endif
    allocate (INDx(FCT%NEQns),stat=IRR)
    if ( IRR/=0 ) then
      write (*,*) '##Error : Not enough memory'
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
      !stop
    endif
    !rmiv
    ndegl = FCT%NDEg*(FCT%NDEg+1)
    ndegl = ndegl/2
    ndeg2 = FCT%NDEg*FCT%NDEg
    !rmiv
    select case (FCT%NDEg)
      case (1)
        call NUFCT(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,INDx,TEMp,FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop,Ir)
      case (2)
        call NUFCT2(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,INDx,TEMp,FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop,Ir)
      case (3)
        call NUFCT3(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,INDx,TEMp,FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop,Ir)
      case (6)
        call NUFCT6(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,INDx,TEMp,FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop,Ir)
      case default
        call NUFCTX(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,INDx,TEMp,FCT%NEQns,FCT%PARent,FCT%NCH,FCT%NSTop,&
             FCT%NDEg,ndegl,Ir)
    endselect

    FCT%STAge = 30
    deallocate (TEMp)
    deallocate (INDx)
  end subroutine NUFCT0

  !======================================================================!
  !> @brief NUFCT performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !======================================================================!
  subroutine NUFCT(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Ir
    integer(kind=kint), intent(out):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Zln(:)
    real(kind=kreal), intent(inout):: Diag(:)
    real(kind=kreal), intent(out):: Temp(:)
    real(kind=kreal), intent(inout):: Dsln(:)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: l
    integer(kind=kint):: ISEm
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5
    real(kind=kreal):: tt

    ISEm = 1
    !
    ! phase I
    !
    call PTIME(t1)
    Diag(1) = 1.0D0/Diag(1)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    do ic = 2, Nstop - 1
      call sum(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx,ISEm)
    enddo
    !
    ! phase II
    !
    call PTIME(t2)
    do ic = Nstop, Neqns
      call SUM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    enddo
    !
    ! phase III
    !
    call PTIME(t3)
    call SUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    call PTIME(t4)
    call SUM3(Neqns-Nstop+1,Dsln,Diag(Nstop:),Indx,Temp)
    call PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
    return
    Ir = 30
  end subroutine NUFCT

  !======================================================================!
  !> @brief NUFCT2  performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !======================================================================!
  subroutine NUFCT2(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Ir
    integer(kind=kint), intent(out):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Zln(4,*)
    real(kind=kreal), intent(inout):: Diag(3,*)
    real(kind=kreal), intent(out):: Temp(4,*)
    real(kind=kreal), intent(inout):: Dsln(4,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: l
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5
    real(kind=kreal):: tt

    !
    ! phase I
    !
    Ir = 0
    call PTIME(t1)
    if ( Nstop>1 ) call INV2(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    do ic = 2, Nstop - 1
      call S2UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    enddo
    !
    ! phase II
    !
    call PTIME(t2)
    do ic = Nstop, Neqns
      call S2UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    enddo
    !
    ! phase III
    !
    call PTIME(t3)
    call S2UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    call PTIME(t4)
    call S2UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    call PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  end subroutine NUFCT2

  !======================================================================!
  !> @brief NUFCT3 performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !======================================================================!
  subroutine NUFCT3(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Ir
    integer(kind=kint), intent(out):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Zln(9,*)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(out):: Temp(:)
    real(kind=kreal), intent(inout):: Dsln(9,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: l
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5
    real(kind=kreal):: tt

    !
    ! phase I
    !
    call PTIME(t1)
    if ( Nstop>1 ) call INV3(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    do ic = 2, Nstop - 1
      call S3UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    enddo
    !
    ! phase II
    !
    call PTIME(t2)
    do ic = Nstop, Neqns
      call S3UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    enddo
    !
    ! phase III
    !
    call PTIME(t3)
    call S3UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    call PTIME(t4)
    call S3UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    call PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  end subroutine NUFCT3

  !======================================================================!
  !> @brief NUFCT6 performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !======================================================================!
  subroutine NUFCT6(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Zln(36,*)
    real(kind=kreal), intent(inout):: Diag(21,*)
    real(kind=kreal), intent(out):: Temp(:)
    real(kind=kreal), intent(inout):: Dsln(36,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: Ir
    integer(kind=kint):: l
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5
    real(kind=kreal):: tt

    !
    ! phase I
    !
    call PTIME(t1)
    if ( Nstop>1 ) call INV6(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    do ic = 2, Nstop - 1
      call S6UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    enddo
    !
    ! phase II
    !
    call PTIME(t2)
    do ic = Nstop, Neqns
      call S6UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    enddo
    !
    ! phase III
    !
    call PTIME(t3)
    call S6UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    call PTIME(t4)
    call S6UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    call PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  end subroutine NUFCT6

  !======================================================================!
  !> @brief NUFCTX performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !======================================================================!
  subroutine NUFCTX(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ndeg,Ndegl,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Parent(:)
    integer(kind=kint), intent(out):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Zln(Ndeg*Ndeg,*)
    real(kind=kreal), intent(inout):: Diag(Ndegl,*)
    real(kind=kreal), intent(out):: Temp(Ndeg*Ndeg,*)
    real(kind=kreal), intent(inout):: Dsln(Ndeg*Ndeg,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: Ir
    integer(kind=kint):: l
    real(kind=kreal):: t1
    real(kind=kreal):: t2
    real(kind=kreal):: t3
    real(kind=kreal):: t4
    real(kind=kreal):: t5
    real(kind=kreal):: tt
    real(kind=kreal):: zz(100)
    real(kind=kreal):: t(100)

    !
    ! phase I
    !
    call PTIME(t1)
    if ( Nstop>1 ) call INVX(Diag(1,1),Ndeg,Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    do ic = 2, Nstop - 1
      call SXUM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx,Ndeg,Ndegl,zz,t)
    enddo
    !
    ! phase II
    !
    call PTIME(t2)
    do ic = Nstop, Neqns
      call SXUM1(ic,Xlnzr,Colno,Zln,Temp,Indx,Ndeg,t)
    enddo
    !
    ! phase III
    !
    call PTIME(t3)
    call SXUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx,Ndeg,Ndegl)
    !
    ! phase IV
    !
    call PTIME(t4)
    call SXUM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp,Ndeg,Ndegl,t)
    call PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  end subroutine NUFCTX

  !======================================================================!
  !> @brief SUM
  !======================================================================!
  subroutine sum(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx,ISEm)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Par(:)
    integer(kind=kint), intent(inout):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    integer(kind=kint), intent(inout):: ISEm
    real(kind=kreal), intent(inout):: Diag(:)
    real(kind=kreal), intent(inout):: Temp(:)
    real(kind=kreal), intent(inout):: Zln(:)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    real(kind=kreal):: piv
    real(kind=kreal):: s
    real(kind=kreal):: t
    real(kind=kreal):: zz

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t = 0.0D0

    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      s = 0.0D0
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) s = s + Temp(j)*Zln(jj)
      enddo
      zz = Zln(k) - s
      Zln(k) = zz*Diag(jc)
      Temp(jc) = zz
      t = t + zz*Zln(k)
    enddo
    piv = Diag(Ic) - t
    if ( dabs(piv)>RMIn ) Diag(Ic) = 1.0D0/piv
    do while ( ISEm/=1 )
    enddo
    ISEm = 0
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
    ISEm = 1
  end subroutine sum

  !======================================================================!
  !> @brief SUM1
  !======================================================================!
  subroutine SUM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Temp(:)
    real(kind=kreal), intent(inout):: Zln(:)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    real(kind=kreal):: s
    real(kind=kreal):: t
    real(kind=kreal):: zz

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t = 0.0D0

    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      s = 0.0D0
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) s = s + Temp(j)*Zln(jj)
      enddo
      zz = Zln(k) - s

      Zln(k) = zz
      Temp(jc) = zz
    enddo
  end subroutine SUM1

  !======================================================================!
  !> @brief SUM2
  !======================================================================!
  subroutine SUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(:)
    real(kind=kreal), intent(inout):: Dsln(:)
    real(kind=kreal), intent(inout):: Temp(:)
    real(kind=kreal), intent(inout):: Zln(:)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    real(kind=kreal):: s

    joc = 0
    do ic = Nstop, Neqns
      Temp(1:Nstop) = 0.0D0
      ks = Xlnzr(ic)
      ke = Xlnzr(ic+1) - 1
      do k = ks, ke
        jj = Colno(k)
        Temp(jj) = Zln(k)
        Zln(k) = Temp(jj)*Diag(jj)
        Indx(jj) = ic
        Diag(ic) = Diag(ic) - Temp(jj)*Zln(k)
      enddo
      do jc = Nstop, ic - 1
        s = 0.0D0
        joc = joc + 1
        do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          if ( Indx(j)==ic ) s = s + Temp(j)*Zln(jj)
        enddo
        if ( s==0.0D0 ) write (16,*) ic, jc
        Dsln(joc) = Dsln(joc) - s
      enddo
    enddo
  end subroutine SUM2

  !======================================================================!
  !> @brief SUM3
  !======================================================================!
  subroutine SUM3(N,Dsln,Diag,Indx,Temp)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    integer(kind=kint), intent(out):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(:)
    real(kind=kreal), intent(inout):: Dsln(:)
    real(kind=kreal), intent(out):: Temp(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc

    if ( N>0 ) then
      Indx(1) = 0
      joc = 1
      Diag(1) = 1.0D0/Diag(1)
      do i = 2, N
        Indx(i) = joc
        do j = 1, i - 1
          Dsln(joc) = Dsln(joc) - DDOT(Dsln(Indx(i):),Dsln(Indx(j):),j-1)
          joc = joc + 1
        enddo
        call VPROD(Dsln(Indx(i):),Diag,Temp,i-1)
        Diag(i) = Diag(i) - DDOT(Temp,Dsln(Indx(i):),i-1)
        call VCOPY(Temp,Dsln(Indx(i):),i-1)
        Diag(i) = 1.0D0/Diag(i)
      enddo
    endif
  end subroutine SUM3

  !======================================================================!
  !> @brief S2UM
  !======================================================================!
  subroutine S2UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Par(:)
    integer(kind=kint), intent(inout):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Diag(3,*)
    real(kind=kreal), intent(inout):: Temp(4,*)
    real(kind=kreal), intent(inout):: Zln(4,*)
    !------
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    real(kind=kreal):: s(4)
    real(kind=kreal):: t(3)
    real(kind=kreal):: zz(4)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t(1:3) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      s(1:4) = 0.0D0
      zz(1:4) = Zln(1:4,k)
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(3,j)*Zln(3,jj)
          zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(4,j)*Zln(3,jj)
          zz(3) = zz(3) - Temp(1,j)*Zln(2,jj) - Temp(3,j)*Zln(4,jj)
          zz(4) = zz(4) - Temp(2,j)*Zln(2,jj) - Temp(4,j)*Zln(4,jj)
        endif
      enddo
      call INV22(Zln(1,k),zz,Diag(1,jc))
      Temp(1:4,jc) = zz(1:4)
      t(1) = t(1) + zz(1)*Zln(1,k) + zz(3)*Zln(3,k)
      t(2) = t(2) + zz(1)*Zln(2,k) + zz(3)*Zln(4,k)
      t(3) = t(3) + zz(2)*Zln(2,k) + zz(4)*Zln(4,k)
    enddo
    Diag(1,Ic) = Diag(1,Ic) - t(1)
    Diag(2,Ic) = Diag(2,Ic) - t(2)
    Diag(3,Ic) = Diag(3,Ic) - t(3)
    call INV2(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  end subroutine S2UM

  !======================================================================!
  !> @brief S2UM1
  !======================================================================!
  subroutine S2UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Temp(4,*)
    real(kind=kreal), intent(inout):: Zln(4,*)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    real(kind=kreal):: s(4)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    s(1:4) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(3,j)*Zln(3,jj)
          s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(4,j)*Zln(3,jj)
          s(3) = s(3) + Temp(1,j)*Zln(2,jj) + Temp(3,j)*Zln(4,jj)
          s(4) = s(4) + Temp(2,j)*Zln(2,jj) + Temp(4,j)*Zln(4,jj)
        endif
      enddo
      do l = 1, 4
        Temp(l,jc) = Zln(l,k) - s(l)
        Zln(l,k) = Temp(l,jc)
        s(l) = 0.0D0
      enddo
    enddo
  end subroutine S2UM1

  !======================================================================!
  !> @brief S2UM2
  !======================================================================!
  subroutine S2UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(3,*)
    real(kind=kreal), intent(inout):: Dsln(4,*)
    real(kind=kreal), intent(inout):: Temp(4,*)
    real(kind=kreal), intent(inout):: Zln(4,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    joc = 0
    do ic = Nstop, Neqns
      ks = Xlnzr(ic)
      ke = Xlnzr(ic+1) - 1
      do k = ks, ke
        jj = Colno(k)
        Temp(1,jj) = Zln(1,k)
        Temp(2,jj) = Zln(2,k)
        Temp(3,jj) = Zln(3,k)
        Temp(4,jj) = Zln(4,k)

        Zln(3,k) = Temp(3,jj) - Temp(1,jj)*Diag(2,jj)
        Zln(1,k) = Temp(1,jj)*Diag(1,jj)
        Zln(3,k) = Zln(3,k)*Diag(3,jj)
        Zln(1,k) = Zln(1,k) - Zln(3,k)*Diag(2,jj)

        Zln(4,k) = Temp(4,jj) - Temp(2,jj)*Diag(2,jj)
        Zln(2,k) = Temp(2,jj)*Diag(1,jj)
        Zln(4,k) = Zln(4,k)*Diag(3,jj)
        Zln(2,k) = Zln(2,k) - Zln(4,k)*Diag(2,jj)

        Diag(1,ic) = Diag(1,ic) - (Temp(1,jj)*Zln(1,k)+Temp(3,jj)*Zln(3,k))
        Diag(2,ic) = Diag(2,ic) - (Temp(1,jj)*Zln(2,k)+Temp(3,jj)*Zln(4,k))
        Diag(3,ic) = Diag(3,ic) - (Temp(2,jj)*Zln(2,k)+Temp(4,jj)*Zln(4,k))
        Indx(jj) = ic
      enddo
      do jc = Nstop, ic - 1
        joc = joc + 1
        do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          if ( Indx(j)==ic ) then
            Dsln(1,joc) = Dsln(1,joc) - (Temp(1,j)*Zln(1,jj)+Temp(3,j)*Zln(3,jj))
            Dsln(2,joc) = Dsln(2,joc) - (Temp(2,j)*Zln(1,jj)+Temp(4,j)*Zln(3,jj))
            Dsln(3,joc) = Dsln(3,joc) - (Temp(1,j)*Zln(2,jj)+Temp(3,j)*Zln(4,jj))
            Dsln(4,joc) = Dsln(4,joc) - (Temp(2,j)*Zln(2,jj)+Temp(4,j)*Zln(4,jj))
          endif
        enddo
      enddo
    enddo
  end subroutine S2UM2

  !======================================================================!
  !> @brief S2UM3
  !======================================================================!
  subroutine S2UM3(N,Dsln,Diag,Indx,Temp)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    integer(kind=kint), intent(out):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(3,*)
    real(kind=kreal), intent(inout):: Dsln(4,*)
    real(kind=kreal), intent(out):: Temp(4,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: joc
    real(kind=kreal):: t(4)

    if ( N>0 ) then
      Indx(1) = 0
      joc = 1
      call INV2(Diag(1,1),ir)
      do i = 2, N
        Indx(i) = joc
        do j = 1, i - 1
          call D2DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
          Dsln(1,joc) = Dsln(1,joc) - t(1)
          Dsln(2,joc) = Dsln(2,joc) - t(2)
          Dsln(3,joc) = Dsln(3,joc) - t(3)
          Dsln(4,joc) = Dsln(4,joc) - t(4)
          joc = joc + 1
        enddo
        call V2PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
        call D2DOT(t,Temp,Dsln(1,Indx(i)),i-1)
        Diag(1,i) = Diag(1,i) - t(1)
        Diag(2,i) = Diag(2,i) - t(2)
        Diag(3,i) = Diag(3,i) - t(4)
        call VCOPY(Temp,Dsln(1,Indx(i)),4*(i-1))
        call INV2(Diag(1,i),ir)
      enddo
    endif
  end subroutine S2UM3

  !======================================================================!
  !> @brief S3UM
  !======================================================================!
  subroutine S3UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Par(:)
    integer(kind=kint), intent(inout):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(inout):: Temp(9,*)
    real(kind=kreal), intent(inout):: Zln(9,*)
    !------
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    integer(kind=kint):: l
    real(kind=kreal):: t(6)
    real(kind=kreal):: zz(9)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t(1:6) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      zz(1:9) = Zln(1:9,k)
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(4,j)*Zln(4,jj) - Temp(7,j)*Zln(7,jj)
          zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(5,j)*Zln(4,jj) - Temp(8,j)*Zln(7,jj)
          zz(3) = zz(3) - Temp(3,j)*Zln(1,jj) - Temp(6,j)*Zln(4,jj) - Temp(9,j)*Zln(7,jj)

          zz(4) = zz(4) - Temp(1,j)*Zln(2,jj) - Temp(4,j)*Zln(5,jj) - Temp(7,j)*Zln(8,jj)
          zz(5) = zz(5) - Temp(2,j)*Zln(2,jj) - Temp(5,j)*Zln(5,jj) - Temp(8,j)*Zln(8,jj)
          zz(6) = zz(6) - Temp(3,j)*Zln(2,jj) - Temp(6,j)*Zln(5,jj) - Temp(9,j)*Zln(8,jj)

          zz(7) = zz(7) - Temp(1,j)*Zln(3,jj) - Temp(4,j)*Zln(6,jj) - Temp(7,j)*Zln(9,jj)
          zz(8) = zz(8) - Temp(2,j)*Zln(3,jj) - Temp(5,j)*Zln(6,jj) - Temp(8,j)*Zln(9,jj)
          zz(9) = zz(9) - Temp(3,j)*Zln(3,jj) - Temp(6,j)*Zln(6,jj) - Temp(9,j)*Zln(9,jj)
        endif
      enddo
      call INV33(Zln(1,k),zz,Diag(1,jc))
      Temp(1:9,jc) = zz(1:9)
      t(1) = t(1) + zz(1)*Zln(1,k) + zz(4)*Zln(4,k) + zz(7)*Zln(7,k)
      t(2) = t(2) + zz(1)*Zln(2,k) + zz(4)*Zln(5,k) + zz(7)*Zln(8,k)
      t(3) = t(3) + zz(2)*Zln(2,k) + zz(5)*Zln(5,k) + zz(8)*Zln(8,k)
      t(4) = t(4) + zz(1)*Zln(3,k) + zz(4)*Zln(6,k) + zz(7)*Zln(9,k)
      t(5) = t(5) + zz(2)*Zln(3,k) + zz(5)*Zln(6,k) + zz(8)*Zln(9,k)
      t(6) = t(6) + zz(3)*Zln(3,k) + zz(6)*Zln(6,k) + zz(9)*Zln(9,k)
    enddo
    do l = 1, 6
      Diag(l,Ic) = Diag(l,Ic) - t(l)
    enddo
    call INV3(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  end subroutine S3UM

  !======================================================================!
  !> @brief S3UM1
  !======================================================================!
  subroutine S3UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Temp(9,*)
    real(kind=kreal), intent(inout):: Zln(9,*)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    real(kind=kreal):: s(9)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    s(1:9) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(4,j)*Zln(4,jj) + Temp(7,j)*Zln(7,jj)
          s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(5,j)*Zln(4,jj) + Temp(8,j)*Zln(7,jj)
          s(3) = s(3) + Temp(3,j)*Zln(1,jj) + Temp(6,j)*Zln(4,jj) + Temp(9,j)*Zln(7,jj)

          s(4) = s(4) + Temp(1,j)*Zln(2,jj) + Temp(4,j)*Zln(5,jj) + Temp(7,j)*Zln(8,jj)
          s(5) = s(5) + Temp(2,j)*Zln(2,jj) + Temp(5,j)*Zln(5,jj) + Temp(8,j)*Zln(8,jj)
          s(6) = s(6) + Temp(3,j)*Zln(2,jj) + Temp(6,j)*Zln(5,jj) + Temp(9,j)*Zln(8,jj)

          s(7) = s(7) + Temp(1,j)*Zln(3,jj) + Temp(4,j)*Zln(6,jj) + Temp(7,j)*Zln(9,jj)
          s(8) = s(8) + Temp(2,j)*Zln(3,jj) + Temp(5,j)*Zln(6,jj) + Temp(8,j)*Zln(9,jj)
          s(9) = s(9) + Temp(3,j)*Zln(3,jj) + Temp(6,j)*Zln(6,jj) + Temp(9,j)*Zln(9,jj)
        endif
      enddo
      do l = 1, 9
        Temp(l,jc) = Zln(l,k) - s(l)
        Zln(l,k) = Temp(l,jc)
        s(l) = 0.0D0
      enddo
    enddo
  end subroutine S3UM1

  !======================================================================!
  !> @brief S3UM1
  !======================================================================!
  subroutine S3UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(inout):: Dsln(9,*)
    real(kind=kreal), intent(inout):: Temp(Neqns,9)
    real(kind=kreal), intent(inout):: Zln(9,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: j
    integer(kind=kint):: j1
    integer(kind=kint):: j2
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    joc = 0
    do ic = Nstop, Neqns
      ks = Xlnzr(ic)
      ke = Xlnzr(ic+1) - 1
      do k = ks, ke
        jj = Colno(k)
        Temp(jj,1) = Zln(1,k)
        Temp(jj,2) = Zln(2,k)
        Temp(jj,3) = Zln(3,k)
        Temp(jj,4) = Zln(4,k)
        Temp(jj,5) = Zln(5,k)
        Temp(jj,6) = Zln(6,k)
        Temp(jj,7) = Zln(7,k)
        Temp(jj,8) = Zln(8,k)
        Temp(jj,9) = Zln(9,k)
        Indx(jj) = ic
      enddo

      do k = ks, ke
        jj = Colno(k)
        Zln(4,k) = Temp(jj,4) - Temp(jj,1)*Diag(2,jj)
        Zln(7,k) = Temp(jj,7) - Temp(jj,1)*Diag(4,jj) - Zln(4,k)*Diag(5,jj)
        Zln(1,k) = Temp(jj,1)*Diag(1,jj)
        Zln(4,k) = Zln(4,k)*Diag(3,jj)
        Zln(7,k) = Zln(7,k)*Diag(6,jj)
        Zln(4,k) = Zln(4,k) - Zln(7,k)*Diag(5,jj)
        Zln(1,k) = Zln(1,k) - Zln(4,k)*Diag(2,jj) - Zln(7,k)*Diag(4,jj)

        Zln(5,k) = Temp(jj,5) - Temp(jj,2)*Diag(2,jj)
        Zln(8,k) = Temp(jj,8) - Temp(jj,2)*Diag(4,jj) - Zln(5,k)*Diag(5,jj)
        Zln(2,k) = Temp(jj,2)*Diag(1,jj)
        Zln(5,k) = Zln(5,k)*Diag(3,jj)
        Zln(8,k) = Zln(8,k)*Diag(6,jj)
        Zln(5,k) = Zln(5,k) - Zln(8,k)*Diag(5,jj)
        Zln(2,k) = Zln(2,k) - Zln(5,k)*Diag(2,jj) - Zln(8,k)*Diag(4,jj)

        Zln(6,k) = Temp(jj,6) - Temp(jj,3)*Diag(2,jj)
        Zln(9,k) = Temp(jj,9) - Temp(jj,3)*Diag(4,jj) - Zln(6,k)*Diag(5,jj)
        Zln(3,k) = Temp(jj,3)*Diag(1,jj)
        Zln(6,k) = Zln(6,k)*Diag(3,jj)
        Zln(9,k) = Zln(9,k)*Diag(6,jj)
        Zln(6,k) = Zln(6,k) - Zln(9,k)*Diag(5,jj)
        Zln(3,k) = Zln(3,k) - Zln(6,k)*Diag(2,jj) - Zln(9,k)*Diag(4,jj)
      enddo

      do k = ks, ke
        jj = Colno(k)
        Diag(1,ic) = Diag(1,ic) - Temp(jj,1)*Zln(1,k) - Temp(jj,4)*Zln(4,k) - Temp(jj,7)*Zln(7,k)
        Diag(2,ic) = Diag(2,ic) - Temp(jj,1)*Zln(2,k) - Temp(jj,4)*Zln(5,k) - Temp(jj,7)*Zln(8,k)
        Diag(3,ic) = Diag(3,ic) - Temp(jj,2)*Zln(2,k) - Temp(jj,5)*Zln(5,k) - Temp(jj,8)*Zln(8,k)
        Diag(4,ic) = Diag(4,ic) - Temp(jj,1)*Zln(3,k) - Temp(jj,4)*Zln(6,k) - Temp(jj,7)*Zln(9,k)
        Diag(5,ic) = Diag(5,ic) - Temp(jj,2)*Zln(3,k) - Temp(jj,5)*Zln(6,k) - Temp(jj,8)*Zln(9,k)
        Diag(6,ic) = Diag(6,ic) - Temp(jj,3)*Zln(3,k) - Temp(jj,6)*Zln(6,k) - Temp(jj,9)*Zln(9,k)
      enddo

      do jc = Nstop, ic - 1
        joc = joc + 1
        j1 = Xlnzr(jc)
        j2 = Xlnzr(jc+1)
        do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          if ( Indx(j)==ic ) then
            Dsln(1,joc) = Dsln(1,joc) - Temp(j,1)*Zln(1,jj) - Temp(j,4)*Zln(4,jj) - Temp(j,7)*Zln(7,jj)
            Dsln(2,joc) = Dsln(2,joc) - Temp(j,2)*Zln(1,jj) - Temp(j,5)*Zln(4,jj) - Temp(j,8)*Zln(7,jj)
            Dsln(3,joc) = Dsln(3,joc) - Temp(j,3)*Zln(1,jj) - Temp(j,6)*Zln(4,jj) - Temp(j,9)*Zln(7,jj)

            Dsln(4,joc) = Dsln(4,joc) - Temp(j,1)*Zln(2,jj) - Temp(j,4)*Zln(5,jj) - Temp(j,7)*Zln(8,jj)
            Dsln(5,joc) = Dsln(5,joc) - Temp(j,2)*Zln(2,jj) - Temp(j,5)*Zln(5,jj) - Temp(j,8)*Zln(8,jj)
            Dsln(6,joc) = Dsln(6,joc) - Temp(j,3)*Zln(2,jj) - Temp(j,6)*Zln(5,jj) - Temp(j,9)*Zln(8,jj)

            Dsln(7,joc) = Dsln(7,joc) - Temp(j,1)*Zln(3,jj) - Temp(j,4)*Zln(6,jj) - Temp(j,7)*Zln(9,jj)
            Dsln(8,joc) = Dsln(8,joc) - Temp(j,2)*Zln(3,jj) - Temp(j,5)*Zln(6,jj) - Temp(j,8)*Zln(9,jj)
            Dsln(9,joc) = Dsln(9,joc) - Temp(j,3)*Zln(3,jj) - Temp(j,6)*Zln(6,jj) - Temp(j,9)*Zln(9,jj)
          endif
        enddo
      enddo
    enddo
  end subroutine S3UM2

  !======================================================================!
  !> @brief S3UM3
  !======================================================================!
  subroutine S3UM3(N,Dsln,Diag,Indx,Temp)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    integer(kind=kint), intent(out):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(inout):: Dsln(9,*)
    real(kind=kreal), intent(out):: Temp(9,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: joc
    real(kind=kreal):: t(9)

    if ( N>0 ) then
      Indx(1) = 0
      joc = 1
      call INV3(Diag(1,1),ir)
      do i = 2, N
        Indx(i) = joc
        do j = 1, i - 1
          call D3DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
          Dsln(:,joc) = Dsln(:,joc) - t(:)
          joc = joc + 1
        enddo
        call V3PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
        call D3DOTL(t,Temp,Dsln(1,Indx(i)),i-1)
        Diag(:,i) = Diag(:,i) - t(1:6)
        call VCOPY(Temp,Dsln(1,Indx(i)),9*(i-1))
        call INV3(Diag(1,i),ir)
      enddo
    endif
  end subroutine S3UM3

  !======================================================================!
  !> @brief S6UM
  !======================================================================!
  subroutine S6UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Par(:)
    integer(kind=kint), intent(inout):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Diag(21,*)
    real(kind=kreal), intent(inout):: Temp(36,*)
    real(kind=kreal), intent(inout):: Zln(36,*)
    !------
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    integer(kind=kint):: l
    real(kind=kreal):: t(21)
    real(kind=kreal):: zz(36)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t(1:21) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      zz(1:36) = Zln(1:36,k)
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(7,j)*Zln(7,jj)&
            - Temp(13,j)*Zln(13,jj) - Temp(19,j)*Zln(19,jj)&
            - Temp(25,j)*Zln(25,jj) - Temp(31,j)*Zln(31,jj)
          zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(8,j)*Zln(7,jj)&
            - Temp(14,j)*Zln(13,jj) - Temp(20,j)*Zln(19,jj)&
            - Temp(26,j)*Zln(25,jj) - Temp(32,j)*Zln(31,jj)
          zz(3) = zz(3) - Temp(3,j)*Zln(1,jj) - Temp(9,j)*Zln(7,jj)&
            - Temp(15,j)*Zln(13,jj) - Temp(21,j)*Zln(19,jj)&
            - Temp(27,j)*Zln(25,jj) - Temp(33,j)*Zln(31,jj)
          zz(4) = zz(4) - Temp(4,j)*Zln(1,jj) - Temp(10,j)*Zln(7,jj)&
            - Temp(16,j)*Zln(13,jj) - Temp(22,j)*Zln(19,jj)&
            - Temp(28,j)*Zln(25,jj) - Temp(34,j)*Zln(31,jj)
          zz(5) = zz(5) - Temp(5,j)*Zln(1,jj) - Temp(11,j)*Zln(7,jj)&
            - Temp(17,j)*Zln(13,jj) - Temp(23,j)*Zln(19,jj)&
            - Temp(29,j)*Zln(25,jj) - Temp(35,j)*Zln(31,jj)
          zz(6) = zz(6) - Temp(6,j)*Zln(1,jj) - Temp(12,j)*Zln(7,jj)&
            - Temp(18,j)*Zln(13,jj) - Temp(24,j)*Zln(19,jj)&
            - Temp(30,j)*Zln(25,jj) - Temp(36,j)*Zln(31,jj)
          zz(7) = zz(7) - Temp(1,j)*Zln(2,jj) - Temp(7,j)*Zln(8,jj)&
            - Temp(13,j)*Zln(14,jj) - Temp(19,j)*Zln(20,jj)&
            - Temp(25,j)*Zln(26,jj) - Temp(31,j)*Zln(32,jj)
          zz(8) = zz(8) - Temp(2,j)*Zln(2,jj) - Temp(8,j)*Zln(8,jj)&
            - Temp(14,j)*Zln(14,jj) - Temp(20,j)*Zln(20,jj)&
            - Temp(26,j)*Zln(26,jj) - Temp(32,j)*Zln(32,jj)
          zz(9) = zz(9) - Temp(3,j)*Zln(2,jj) - Temp(9,j)*Zln(8,jj)&
            - Temp(15,j)*Zln(14,jj) - Temp(21,j)*Zln(20,jj)&
            - Temp(27,j)*Zln(26,jj) - Temp(33,j)*Zln(32,jj)
          zz(10) = zz(10) - Temp(4,j)*Zln(2,jj) - Temp(10,j)*Zln(8,jj)&
            - Temp(16,j)*Zln(14,jj) - Temp(22,j)*Zln(20,jj)&
            - Temp(28,j)*Zln(26,jj) - Temp(34,j)*Zln(32,jj)
          zz(11) = zz(11) - Temp(5,j)*Zln(2,jj) - Temp(11,j)*Zln(8,jj)&
            - Temp(17,j)*Zln(14,jj) - Temp(23,j)*Zln(20,jj)&
            - Temp(29,j)*Zln(26,jj) - Temp(35,j)*Zln(32,jj)
          zz(12) = zz(12) - Temp(6,j)*Zln(2,jj) - Temp(12,j)*Zln(8,jj)&
            - Temp(18,j)*Zln(14,jj) - Temp(24,j)*Zln(20,jj)&
            - Temp(30,j)*Zln(26,jj) - Temp(36,j)*Zln(32,jj)
          zz(13) = zz(13) - Temp(1,j)*Zln(3,jj) - Temp(7,j)*Zln(9,jj)&
            - Temp(13,j)*Zln(15,jj) - Temp(19,j)*Zln(21,jj)&
            - Temp(25,j)*Zln(27,jj) - Temp(31,j)*Zln(33,jj)
          zz(14) = zz(14) - Temp(2,j)*Zln(3,jj) - Temp(8,j)*Zln(9,jj)&
            - Temp(14,j)*Zln(15,jj) - Temp(20,j)*Zln(21,jj)&
            - Temp(26,j)*Zln(27,jj) - Temp(32,j)*Zln(33,jj)
          zz(15) = zz(15) - Temp(3,j)*Zln(3,jj) - Temp(9,j)*Zln(9,jj)&
            - Temp(15,j)*Zln(15,jj) - Temp(21,j)*Zln(21,jj)&
            - Temp(27,j)*Zln(27,jj) - Temp(33,j)*Zln(33,jj)
          zz(16) = zz(16) - Temp(4,j)*Zln(3,jj) - Temp(10,j)*Zln(9,jj)&
            - Temp(16,j)*Zln(15,jj) - Temp(22,j)*Zln(21,jj)&
            - Temp(28,j)*Zln(27,jj) - Temp(34,j)*Zln(33,jj)
          zz(17) = zz(17) - Temp(5,j)*Zln(3,jj) - Temp(11,j)*Zln(9,jj)&
            - Temp(17,j)*Zln(15,jj) - Temp(23,j)*Zln(21,jj)&
            - Temp(29,j)*Zln(27,jj) - Temp(35,j)*Zln(33,jj)
          zz(18) = zz(18) - Temp(6,j)*Zln(3,jj) - Temp(12,j)*Zln(9,jj)&
            - Temp(18,j)*Zln(15,jj) - Temp(24,j)*Zln(21,jj)&
            - Temp(30,j)*Zln(27,jj) - Temp(36,j)*Zln(33,jj)
          zz(19) = zz(19) - Temp(1,j)*Zln(4,jj) - Temp(7,j)*Zln(10,jj)&
            - Temp(13,j)*Zln(16,jj) - Temp(19,j)*Zln(22,jj)&
            - Temp(25,j)*Zln(28,jj) - Temp(31,j)*Zln(34,jj)
          zz(20) = zz(20) - Temp(2,j)*Zln(4,jj) - Temp(8,j)*Zln(10,jj)&
            - Temp(14,j)*Zln(16,jj) - Temp(20,j)*Zln(22,jj)&
            - Temp(26,j)*Zln(28,jj) - Temp(32,j)*Zln(34,jj)
          zz(21) = zz(21) - Temp(3,j)*Zln(4,jj) - Temp(9,j)*Zln(10,jj)&
            - Temp(15,j)*Zln(16,jj) - Temp(21,j)*Zln(22,jj)&
            - Temp(27,j)*Zln(28,jj) - Temp(33,j)*Zln(34,jj)
          zz(22) = zz(22) - Temp(4,j)*Zln(4,jj) - Temp(10,j)*Zln(10,jj)&
            - Temp(16,j)*Zln(16,jj) - Temp(22,j)*Zln(22,jj)&
            - Temp(28,j)*Zln(28,jj) - Temp(34,j)*Zln(34,jj)
          zz(23) = zz(23) - Temp(5,j)*Zln(4,jj) - Temp(11,j)*Zln(10,jj)&
            - Temp(17,j)*Zln(16,jj) - Temp(23,j)*Zln(22,jj)&
            - Temp(29,j)*Zln(28,jj) - Temp(35,j)*Zln(34,jj)
          zz(24) = zz(24) - Temp(6,j)*Zln(4,jj) - Temp(12,j)*Zln(10,jj)&
            - Temp(18,j)*Zln(16,jj) - Temp(24,j)*Zln(22,jj)&
            - Temp(30,j)*Zln(28,jj) - Temp(36,j)*Zln(34,jj)
          zz(25) = zz(25) - Temp(1,j)*Zln(5,jj) - Temp(7,j)*Zln(11,jj)&
            - Temp(13,j)*Zln(17,jj) - Temp(19,j)*Zln(23,jj)&
            - Temp(25,j)*Zln(29,jj) - Temp(31,j)*Zln(35,jj)
          zz(26) = zz(26) - Temp(2,j)*Zln(5,jj) - Temp(8,j)*Zln(11,jj)&
            - Temp(14,j)*Zln(17,jj) - Temp(20,j)*Zln(23,jj)&
            - Temp(26,j)*Zln(29,jj) - Temp(32,j)*Zln(35,jj)
          zz(27) = zz(27) - Temp(3,j)*Zln(5,jj) - Temp(9,j)*Zln(11,jj)&
            - Temp(15,j)*Zln(17,jj) - Temp(21,j)*Zln(23,jj)&
            - Temp(27,j)*Zln(29,jj) - Temp(33,j)*Zln(35,jj)
          zz(28) = zz(28) - Temp(4,j)*Zln(5,jj) - Temp(10,j)*Zln(11,jj)&
            - Temp(16,j)*Zln(17,jj) - Temp(22,j)*Zln(23,jj)&
            - Temp(28,j)*Zln(29,jj) - Temp(34,j)*Zln(35,jj)
          zz(29) = zz(29) - Temp(5,j)*Zln(5,jj) - Temp(11,j)*Zln(11,jj)&
            - Temp(17,j)*Zln(17,jj) - Temp(23,j)*Zln(23,jj)&
            - Temp(29,j)*Zln(29,jj) - Temp(35,j)*Zln(35,jj)
          zz(30) = zz(30) - Temp(6,j)*Zln(5,jj) - Temp(12,j)*Zln(11,jj)&
            - Temp(18,j)*Zln(17,jj) - Temp(24,j)*Zln(23,jj)&
            - Temp(30,j)*Zln(29,jj) - Temp(36,j)*Zln(35,jj)
          zz(31) = zz(31) - Temp(1,j)*Zln(6,jj) - Temp(7,j)*Zln(12,jj)&
            - Temp(13,j)*Zln(18,jj) - Temp(19,j)*Zln(24,jj)&
            - Temp(25,j)*Zln(30,jj) - Temp(31,j)*Zln(36,jj)
          zz(32) = zz(32) - Temp(2,j)*Zln(6,jj) - Temp(8,j)*Zln(12,jj)&
            - Temp(14,j)*Zln(18,jj) - Temp(20,j)*Zln(24,jj)&
            - Temp(26,j)*Zln(30,jj) - Temp(32,j)*Zln(36,jj)
          zz(33) = zz(33) - Temp(3,j)*Zln(6,jj) - Temp(9,j)*Zln(12,jj)&
            - Temp(15,j)*Zln(18,jj) - Temp(21,j)*Zln(24,jj)&
            - Temp(27,j)*Zln(30,jj) - Temp(33,j)*Zln(36,jj)
          zz(34) = zz(34) - Temp(4,j)*Zln(6,jj) - Temp(10,j)*Zln(12,jj)&
            - Temp(16,j)*Zln(18,jj) - Temp(22,j)*Zln(24,jj)&
            - Temp(28,j)*Zln(30,jj) - Temp(34,j)*Zln(36,jj)
          zz(35) = zz(35) - Temp(5,j)*Zln(6,jj) - Temp(11,j)*Zln(12,jj)&
            - Temp(17,j)*Zln(18,jj) - Temp(23,j)*Zln(24,jj)&
            - Temp(29,j)*Zln(30,jj) - Temp(35,j)*Zln(36,jj)
          zz(36) = zz(36) - Temp(6,j)*Zln(6,jj) - Temp(12,j)*Zln(12,jj)&
            - Temp(18,j)*Zln(18,jj) - Temp(24,j)*Zln(24,jj)&
            - Temp(30,j)*Zln(30,jj) - Temp(36,j)*Zln(36,jj)
        endif
      enddo
      call INV66(Zln(1,k),zz,Diag(1,jc))
      Temp(1:36,jc) = zz(1:36)

      t(1) = t(1) + zz(1)*Zln(1,k) + zz(7)*Zln(7,k) + zz(13)*Zln(13,k)&
        + zz(19)*Zln(19,k) + zz(25)*Zln(25,k) + zz(31)*Zln(31,k)
      t(2) = t(2) + zz(1)*Zln(2,k) + zz(7)*Zln(8,k) + zz(13)*Zln(14,k)&
        + zz(19)*Zln(20,k) + zz(25)*Zln(26,k) + zz(31)*Zln(32,k)
      t(3) = t(3) + zz(2)*Zln(2,k) + zz(8)*Zln(8,k) + zz(14)*Zln(14,k)&
        + zz(20)*Zln(20,k) + zz(26)*Zln(26,k) + zz(32)*Zln(32,k)
      t(4) = t(4) + zz(1)*Zln(3,k) + zz(7)*Zln(9,k) + zz(13)*Zln(15,k)&
        + zz(19)*Zln(21,k) + zz(25)*Zln(27,k) + zz(31)*Zln(33,k)
      t(5) = t(5) + zz(2)*Zln(3,k) + zz(8)*Zln(9,k) + zz(14)*Zln(15,k)&
        + zz(20)*Zln(21,k) + zz(26)*Zln(27,k) + zz(32)*Zln(33,k)
      t(6) = t(6) + zz(3)*Zln(3,k) + zz(9)*Zln(9,k) + zz(15)*Zln(15,k)&
        + zz(21)*Zln(21,k) + zz(27)*Zln(27,k) + zz(33)*Zln(33,k)
      t(7) = t(7) + zz(1)*Zln(4,k) + zz(7)*Zln(10,k) + zz(13)*Zln(16,k)&
        + zz(19)*Zln(22,k) + zz(25)*Zln(28,k) + zz(31)*Zln(34,k)
      t(8) = t(8) + zz(2)*Zln(4,k) + zz(8)*Zln(10,k) + zz(14)*Zln(16,k)&
        + zz(20)*Zln(22,k) + zz(26)*Zln(28,k) + zz(32)*Zln(34,k)
      t(9) = t(9) + zz(3)*Zln(4,k) + zz(9)*Zln(10,k) + zz(15)*Zln(16,k)&
        + zz(21)*Zln(22,k) + zz(27)*Zln(28,k) + zz(33)*Zln(34,k)
      t(10) = t(10) + zz(4)*Zln(4,k) + zz(10)*Zln(10,k) + zz(16)*Zln(16,k)&
        + zz(22)*Zln(22,k) + zz(28)*Zln(28,k) + zz(34)*Zln(34,k)
      t(11) = t(11) + zz(1)*Zln(5,k) + zz(7)*Zln(11,k) + zz(13)*Zln(17,k)&
        + zz(19)*Zln(23,k) + zz(25)*Zln(29,k) + zz(31)*Zln(35,k)
      t(12) = t(12) + zz(2)*Zln(5,k) + zz(8)*Zln(11,k) + zz(14)*Zln(17,k)&
        + zz(20)*Zln(23,k) + zz(26)*Zln(29,k) + zz(32)*Zln(35,k)
      t(13) = t(13) + zz(3)*Zln(5,k) + zz(9)*Zln(11,k) + zz(15)*Zln(17,k)&
        + zz(21)*Zln(23,k) + zz(27)*Zln(29,k) + zz(33)*Zln(35,k)
      t(14) = t(14) + zz(4)*Zln(5,k) + zz(10)*Zln(11,k) + zz(16)*Zln(17,k)&
        + zz(22)*Zln(23,k) + zz(28)*Zln(29,k) + zz(34)*Zln(35,k)
      t(15) = t(15) + zz(5)*Zln(5,k) + zz(11)*Zln(11,k) + zz(17)*Zln(17,k)&
        + zz(23)*Zln(23,k) + zz(29)*Zln(29,k) + zz(35)*Zln(35,k)
      t(16) = t(16) + zz(1)*Zln(6,k) + zz(7)*Zln(12,k) + zz(13)*Zln(18,k)&
        + zz(19)*Zln(24,k) + zz(25)*Zln(30,k) + zz(31)*Zln(36,k)
      t(17) = t(17) + zz(2)*Zln(6,k) + zz(8)*Zln(12,k) + zz(14)*Zln(18,k)&
        + zz(20)*Zln(24,k) + zz(26)*Zln(30,k) + zz(32)*Zln(36,k)
      t(18) = t(18) + zz(3)*Zln(6,k) + zz(9)*Zln(12,k) + zz(15)*Zln(18,k)&
        + zz(21)*Zln(24,k) + zz(27)*Zln(30,k) + zz(33)*Zln(36,k)
      t(19) = t(19) + zz(4)*Zln(6,k) + zz(10)*Zln(12,k) + zz(16)*Zln(18,k)&
        + zz(22)*Zln(24,k) + zz(28)*Zln(30,k) + zz(34)*Zln(36,k)
      t(20) = t(20) + zz(5)*Zln(6,k) + zz(11)*Zln(12,k) + zz(17)*Zln(18,k)&
        + zz(23)*Zln(24,k) + zz(29)*Zln(30,k) + zz(35)*Zln(36,k)
      t(21) = t(21) + zz(6)*Zln(6,k) + zz(12)*Zln(12,k) + zz(18)*Zln(18,k)&
        + zz(24)*Zln(24,k) + zz(30)*Zln(30,k) + zz(36)*Zln(36,k)
    enddo
    do l = 1, 21
      Diag(l,Ic) = Diag(l,Ic) - t(l)
    enddo
    call INV6(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  end subroutine S6UM

  !======================================================================!
  !> @brief S6UM1
  !======================================================================!
  subroutine S6UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Temp(9,*)
    real(kind=kreal), intent(inout):: Zln(9,*)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    real(kind=kreal):: s(9)

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    s(1:9) = 0.0D0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(4,j)*Zln(4,jj) + Temp(7,j)*Zln(7,jj)
          s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(5,j)*Zln(4,jj) + Temp(8,j)*Zln(7,jj)
          s(3) = s(3) + Temp(3,j)*Zln(1,jj) + Temp(6,j)*Zln(4,jj) + Temp(9,j)*Zln(7,jj)

          s(4) = s(4) + Temp(1,j)*Zln(2,jj) + Temp(4,j)*Zln(5,jj) + Temp(7,j)*Zln(8,jj)
          s(5) = s(5) + Temp(2,j)*Zln(2,jj) + Temp(5,j)*Zln(5,jj) + Temp(8,j)*Zln(8,jj)
          s(6) = s(6) + Temp(3,j)*Zln(2,jj) + Temp(6,j)*Zln(5,jj) + Temp(9,j)*Zln(8,jj)

          s(7) = s(7) + Temp(1,j)*Zln(3,jj) + Temp(4,j)*Zln(6,jj) + Temp(7,j)*Zln(9,jj)
          s(8) = s(8) + Temp(2,j)*Zln(3,jj) + Temp(5,j)*Zln(6,jj) + Temp(8,j)*Zln(9,jj)
          s(9) = s(9) + Temp(3,j)*Zln(3,jj) + Temp(6,j)*Zln(6,jj) + Temp(9,j)*Zln(9,jj)
        endif
      enddo
      do l = 1, 9
        Temp(l,jc) = Zln(l,k) - s(l)
        Zln(l,k) = Temp(l,jc)
        s(l) = 0.0D0
      enddo
    enddo
  end subroutine S6UM1

  !======================================================================!
  !> @brief S6UM2
  !======================================================================!
  subroutine S6UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(21,*)
    real(kind=kreal), intent(inout):: Dsln(36,*)
    real(kind=kreal), intent(inout):: Temp(36,Neqns)
    real(kind=kreal), intent(inout):: Zln(36,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: j
    integer(kind=kint):: j1
    integer(kind=kint):: j2
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    joc = 0
    do ic = Nstop, Neqns
      Temp(1:Nstop,1:36) = 0.0D0
      ks = Xlnzr(ic)
      ke = Xlnzr(ic+1) - 1
      do k = ks, ke
        jj = Colno(k)
        Temp(:,jj) = Zln(:,k)
        Indx(jj) = ic
      enddo
      do k = ks, ke
        jj = Colno(k)
        call INV66(Zln(1,k),Temp,Diag(1,jj))
      enddo

      do k = ks, ke
        jj = Colno(k)
        Diag(1,ic) = Diag(1,ic) - Temp(jj,1)*Zln(1,k) - Temp(jj,4)*Zln(4,k) - Temp(jj,7)*Zln(7,k)
        Diag(2,ic) = Diag(2,ic) - Temp(jj,1)*Zln(2,k) - Temp(jj,4)*Zln(5,k) - Temp(jj,7)*Zln(8,k)
        Diag(3,ic) = Diag(3,ic) - Temp(jj,2)*Zln(2,k) - Temp(jj,5)*Zln(5,k) - Temp(jj,8)*Zln(8,k)
        Diag(4,ic) = Diag(4,ic) - Temp(jj,1)*Zln(3,k) - Temp(jj,4)*Zln(6,k) - Temp(jj,7)*Zln(9,k)
        Diag(5,ic) = Diag(5,ic) - Temp(jj,2)*Zln(3,k) - Temp(jj,5)*Zln(6,k) - Temp(jj,8)*Zln(9,k)
        Diag(6,ic) = Diag(6,ic) - Temp(jj,3)*Zln(3,k) - Temp(jj,6)*Zln(6,k) - Temp(jj,9)*Zln(9,k)
      enddo
      do jc = Nstop, ic - 1
        joc = joc + 1
        j1 = Xlnzr(jc)
        j2 = Xlnzr(jc+1)
        do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          if ( Indx(j)==ic ) then
            Dsln(1,joc) = Dsln(1,joc) - Temp(j,1)*Zln(1,jj) - Temp(j,4)*Zln(4,jj) - Temp(j,7)*Zln(7,jj)
            Dsln(2,joc) = Dsln(2,joc) - Temp(j,2)*Zln(1,jj) - Temp(j,5)*Zln(4,jj) - Temp(j,8)*Zln(7,jj)
            Dsln(3,joc) = Dsln(3,joc) - Temp(j,3)*Zln(1,jj) - Temp(j,6)*Zln(4,jj) - Temp(j,9)*Zln(7,jj)

            Dsln(4,joc) = Dsln(4,joc) - Temp(j,1)*Zln(2,jj) - Temp(j,4)*Zln(5,jj) - Temp(j,7)*Zln(8,jj)
            Dsln(5,joc) = Dsln(5,joc) - Temp(j,2)*Zln(2,jj) - Temp(j,5)*Zln(5,jj) - Temp(j,8)*Zln(8,jj)
            Dsln(6,joc) = Dsln(6,joc) - Temp(j,3)*Zln(2,jj) - Temp(j,6)*Zln(5,jj) - Temp(j,9)*Zln(8,jj)

            Dsln(7,joc) = Dsln(7,joc) - Temp(j,1)*Zln(3,jj) - Temp(j,4)*Zln(6,jj) - Temp(j,7)*Zln(9,jj)
            Dsln(8,joc) = Dsln(8,joc) - Temp(j,2)*Zln(3,jj) - Temp(j,5)*Zln(6,jj) - Temp(j,8)*Zln(9,jj)
            Dsln(9,joc) = Dsln(9,joc) - Temp(j,3)*Zln(3,jj) - Temp(j,6)*Zln(6,jj) - Temp(j,9)*Zln(9,jj)
          endif
        enddo
      enddo
    enddo
  end subroutine S6UM2

  !======================================================================!
  !> @brief S6UM3
  !======================================================================!
  subroutine S6UM3(N,Dsln,Diag,Indx,Temp)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    integer(kind=kint), intent(out):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(6,*)
    real(kind=kreal), intent(inout):: Dsln(9,*)
    real(kind=kreal), intent(out):: Temp(9,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: l
    real(kind=kreal):: t(9)

    if ( N>0 ) then
      Indx(1) = 0
      joc = 1
      call INV3(Diag(1,1),ir)
      do i = 2, N
        Indx(i) = joc
        do j = 1, i - 1
          call D3DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
          do l = 1, 9
            Dsln(l,joc) = Dsln(l,joc) - t(l)
          enddo
          joc = joc + 1
        enddo
        call V3PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
        call D3DOTL(t,Temp,Dsln(1,Indx(i)),i-1)
        do l = 1, 6
          Diag(l,i) = Diag(l,i) - t(l)
        enddo
        call VCOPY(Temp,Dsln(1,Indx(i)),9*(i-1))
        call INV3(Diag(1,i),ir)
      enddo
    endif
  end subroutine S6UM3

  !======================================================================!
  !> @brief SXUM
  !======================================================================!
  subroutine SXUM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx,Ndeg,Ndegl,Zz,T)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ic
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Par(:)
    integer(kind=kint), intent(inout):: Indx(:)
    integer(kind=kint), intent(inout):: Nch(:)
    real(kind=kreal), intent(inout):: Diag(Ndegl,*)
    real(kind=kreal), intent(out):: T(Ndegl)
    real(kind=kreal), intent(inout):: Temp(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Zln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(out):: Zz(Ndeg,Ndeg)
    !------
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    integer(kind=kint):: m
    integer(kind=kint):: n
    integer(kind=kint):: ndeg22

    ndeg22 = Ndeg*Ndeg
    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    T = 0.0
    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      Zz = Zln(:,:,k)
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          do m = 1, Ndeg, 2
            do n = 1, Ndeg, 2
              do kk = 1, Ndeg, 2
                Zz(n,m) = Zz(n,m) - Temp(n,kk,j)*Zln(m,kk,jj) - Temp(n,kk+1,j)*Zln(m,kk+1,jj)
                Zz(n,m+1) = Zz(n,m+1) - Temp(n,kk,j)*Zln(m+1,kk,jj) - Temp(n,kk+1,j)*Zln(m+1,kk+1,jj)
                Zz(n+1,m) = Zz(n+1,m) - Temp(n+1,kk,j)*Zln(m,kk,jj) - Temp(n+1,kk+1,j)*Zln(m,kk+1,jj)
                Zz(n+1,m+1) = Zz(n+1,m+1) - Temp(n+1,kk,j)*Zln(m+1,kk,jj) - Temp(n+1,kk+1,j)*Zln(m+1,kk+1,jj)
              enddo
            enddo
          enddo
        endif
      enddo
      call INVXX(Zln(1,1,k),Zz,Diag(1,jc),Ndeg)

      Temp(:,:,jc) = Zz
      joc = 0
      do n = 1, Ndeg
        do m = 1, n
          joc = joc + 1
          do kk = 1, Ndeg, 2
            T(joc) = T(joc) + Zz(n,kk)*Zln(m,kk,k) + Zz(n,kk+1)*Zln(m,kk+1,k)
          enddo
        enddo
      enddo
    enddo

    Diag(:,Ic) = Diag(:,Ic) - T
    call INVX(Diag(1,Ic),Ndeg,ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  end subroutine SXUM

  !======================================================================!
  !> @brief SXUM1
  !======================================================================!
  subroutine SXUM1(Ic,Xlnzr,Colno,Zln,Temp,Indx,Ndeg,S)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: S(Ndeg,Ndeg)
    real(kind=kreal), intent(inout):: Temp(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Zln(Ndeg,Ndeg,*)
    !------
    integer(kind=kint):: Ic
    integer(kind=kint):: j
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    integer(kind=kint):: m
    integer(kind=kint):: n

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    S(1:Ndeg,1:Ndeg) = 0.0D0

    do k = ks, ke - 1
      jc = Colno(k)
      Indx(jc) = Ic
      do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
        j = Colno(jj)
        if ( Indx(j)==Ic ) then
          do m = 1, Ndeg
            do n = 1, Ndeg
              do kk = 1, Ndeg
                S(n,m) = S(n,m) + Temp(n,kk,j)*Zln(m,kk,jj)
              enddo
            enddo
          enddo
        endif
      enddo
      do m = 1, Ndeg
        do n = 1, Ndeg
          Temp(n,m,jc) = Zln(n,m,k) - S(n,m)
          Zln(n,m,k) = Temp(n,m,jc)
          S(n,m) = 0.0D0
        enddo
      enddo
    enddo
  end subroutine SXUM1

  !======================================================================!
  !> @brief SXUM2
  !======================================================================!
  subroutine SXUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx,Ndeg,Ndegl)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(inout):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(Ndegl,*)
    real(kind=kreal), intent(inout):: Dsln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Temp(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Zln(Ndeg,Ndeg,*)
    !------
    integer(kind=kint):: ic
    integer(kind=kint):: j
    integer(kind=kint):: j1
    integer(kind=kint):: j2
    integer(kind=kint):: jc
    integer(kind=kint):: jj
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: kk
    integer(kind=kint):: ks
    integer(kind=kint):: locd
    integer(kind=kint):: m
    integer(kind=kint):: n

    joc = 0
    do ic = Nstop, Neqns
      ks = Xlnzr(ic)
      ke = Xlnzr(ic+1) - 1
      do k = ks, ke
        jj = Colno(k)
        do m = 1, Ndeg
          Temp(1:Ndeg,m,jj) = Zln(1:Ndeg,m,k)
          Indx(jj) = ic
        enddo
      enddo
      do k = ks, ke
        jj = Colno(k)
        call INVXX(Zln(1,1,k),Temp(1,1,jj),Diag(1,jj),Ndeg)
      enddo

      locd = 0
      do n = 1, Ndeg
        do m = 1, n
          locd = locd + 1
          do k = ks, ke
            jj = Colno(k)
            do kk = 1, Ndeg
              Diag(locd,ic) = Diag(locd,ic) - Temp(n,kk,jj)*Zln(m,kk,k)
            enddo
          enddo
        enddo
      enddo
      do jc = Nstop, ic - 1
        joc = joc + 1
        j1 = Xlnzr(jc)
        j2 = Xlnzr(jc+1)
        do jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          if ( Indx(j)==ic ) then
            do m = 1, Ndeg
              do n = 1, Ndeg
                do k = 1, Ndeg
                  Dsln(n,m,joc) = Dsln(n,m,joc) - Temp(n,k,j)*Zln(m,k,jj)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  end subroutine SXUM2

  !======================================================================!
  !> @brief SXUM3
  !======================================================================!
  subroutine SXUM3(Nn,Dsln,Diag,Indx,Temp,Ndeg,Ndegl,T)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    integer(kind=kint), intent(in):: Nn
    integer(kind=kint), intent(out):: Indx(:)
    real(kind=kreal), intent(inout):: Diag(Ndegl,*)
    real(kind=kreal), intent(inout):: Dsln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(out):: T(Ndeg,Ndeg)
    real(kind=kreal), intent(out):: Temp(Ndeg,Ndeg,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: ir
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: locd
    integer(kind=kint):: m
    integer(kind=kint):: n

    if ( Nn>0 ) then
      Indx(1) = 0
      joc = 1
      call INVX(Diag(1,1),Ndeg,ir)
      do i = 2, Nn
        Indx(i) = joc
        do j = 1, i - 1
          call DXDOT(Ndeg,T,Dsln(1,1,Indx(i)),Dsln(1,1,Indx(j)),j-1)
          do m = 1, Ndeg
            Dsln(1:Ndeg,m,joc) = Dsln(1:Ndeg,m,joc) - T(1:Ndeg,m)
          enddo
          joc = joc + 1
        enddo
        call VXPROD(Ndeg,Ndegl,Dsln(1,1,Indx(i)),Diag,Temp,i-1)
        call DXDOTL(Ndeg,T,Temp,Dsln(1,1,Indx(i)),i-1)
        locd = 0
        do n = 1, Ndeg
          do m = 1, n
            locd = locd + 1
            Diag(locd,i) = Diag(locd,i) - T(n,m)
          enddo
        enddo
        call VCOPY(Temp,Dsln(1,1,Indx(i)),Ndeg*Ndeg*(i-1))
        call INVX(Diag(1,i),Ndeg,ir)
      enddo
    endif
  end subroutine SXUM3

  !======================================================================!
  !> @brief INV2
  !======================================================================!
  subroutine INV2(Dsln,Ir)
    implicit none
    !------
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Dsln(3)
    !------
    real(kind=kreal):: t

    Ir = 0
    if ( dabs(Dsln(1))<RMIn ) then
      Ir = 10
      return
    endif
    Dsln(1) = 1.0D0/Dsln(1)
    t = Dsln(2)*Dsln(1)
    Dsln(3) = Dsln(3) - t*Dsln(2)
    Dsln(2) = t
    if ( dabs(Dsln(3))<RMIn ) then
      Ir = 10
      return
    endif
    Dsln(3) = 1.0D0/Dsln(3)
  end subroutine INV2

  !======================================================================!
  !> @brief INV22
  !======================================================================!
  subroutine INV22(Zln,Zz,Diag)
    implicit none
    !------
    real(kind=kreal), intent(in):: Diag(3)
    real(kind=kreal), intent(in):: Zz(4)
    real(kind=kreal), intent(out):: Zln(4)
    !------
    Zln(3) = Zz(3) - Zz(1)*Diag(2)
    Zln(1) = Zz(1)*Diag(1)
    Zln(3) = Zln(3)*Diag(3)
    Zln(1) = Zln(1) - Zln(3)*Diag(2)

    Zln(4) = Zz(4) - Zz(2)*Diag(2)
    Zln(2) = Zz(2)*Diag(1)
    Zln(4) = Zln(4)*Diag(3)
    Zln(2) = Zln(2) - Zln(4)*Diag(2)
  end subroutine INV22

  !======================================================================!
  !> @brief INV3
  !======================================================================!
  subroutine INV3(Dsln,Ir)
    implicit none
    !------
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Dsln(6)
    !------
    real(kind=kreal):: t(2)

    Ir = 0
    do
      if ( dabs(Dsln(1))<RMIn ) exit
      Dsln(1) = 1.0D0/Dsln(1)
      t(1) = Dsln(2)*Dsln(1)
      Dsln(3) = Dsln(3) - t(1)*Dsln(2)
      Dsln(2) = t(1)
      if ( dabs(Dsln(3))<RMIn ) exit
      Dsln(3) = 1.0D0/Dsln(3)
      t(1) = Dsln(4)*Dsln(1)
      Dsln(5) = Dsln(5) - Dsln(2)*Dsln(4)
      t(2) = Dsln(5)*Dsln(3)
      Dsln(6) = Dsln(6) - t(1)*Dsln(4) - t(2)*Dsln(5)
      Dsln(4) = t(1)
      Dsln(5) = t(2)
      if ( dabs(Dsln(6))<RMIn ) exit
      Dsln(6) = 1.0D0/Dsln(6)
      return
    enddo

    Dsln(1) = 1.0D0
    Dsln(2) = 0.0D0
    Dsln(3) = 1.0D0
    Dsln(4) = 0.0D0
    Dsln(5) = 0.0D0
    Dsln(6) = 1.0D0
  end subroutine INV3

  !======================================================================!
  !> @brief INV33
  !======================================================================!
  subroutine INV33(Zln,Zz,Diag)
    implicit none
    !------
    real(kind=kreal), intent(in):: Diag(6)
    real(kind=kreal), intent(in):: Zz(9)
    real(kind=kreal), intent(out):: Zln(9)
    !------
    Zln(4) = Zz(4) - Zz(1)*Diag(2)
    Zln(7) = Zz(7) - Zz(1)*Diag(4) - Zln(4)*Diag(5)
    Zln(1) = Zz(1)*Diag(1)
    Zln(4) = Zln(4)*Diag(3)
    Zln(7) = Zln(7)*Diag(6)
    Zln(4) = Zln(4) - Zln(7)*Diag(5)
    Zln(1) = Zln(1) - Zln(4)*Diag(2) - Zln(7)*Diag(4)

    Zln(5) = Zz(5) - Zz(2)*Diag(2)
    Zln(8) = Zz(8) - Zz(2)*Diag(4) - Zln(5)*Diag(5)
    Zln(2) = Zz(2)*Diag(1)
    Zln(5) = Zln(5)*Diag(3)
    Zln(8) = Zln(8)*Diag(6)
    Zln(5) = Zln(5) - Zln(8)*Diag(5)
    Zln(2) = Zln(2) - Zln(5)*Diag(2) - Zln(8)*Diag(4)

    Zln(6) = Zz(6) - Zz(3)*Diag(2)
    Zln(9) = Zz(9) - Zz(3)*Diag(4) - Zln(6)*Diag(5)
    Zln(3) = Zz(3)*Diag(1)
    Zln(6) = Zln(6)*Diag(3)
    Zln(9) = Zln(9)*Diag(6)
    Zln(6) = Zln(6) - Zln(9)*Diag(5)
    Zln(3) = Zln(3) - Zln(6)*Diag(2) - Zln(9)*Diag(4)
  end subroutine INV33

  !======================================================================!
  !> @brief INV6
  !======================================================================!
  subroutine INV6(Dsln,Ir)
    implicit none
    !------
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Dsln(21)
    !------
    real(kind=kreal):: t(5)

    Ir = 0
    Dsln(1) = 1.0D0/Dsln(1)
    t(1) = Dsln(2)*Dsln(1)
    Dsln(3) = 1.0D0/(Dsln(3)-t(1)*Dsln(2))
    Dsln(2) = t(1)
    Dsln(5) = Dsln(5) - Dsln(4)*Dsln(2)
    t(1) = Dsln(4)*Dsln(1)
    t(2) = Dsln(5)*Dsln(3)
    Dsln(6) = 1.0D0/(Dsln(6)-t(1)*Dsln(4)-t(2)*Dsln(5))
    Dsln(4) = t(1)
    Dsln(5) = t(2)
    Dsln(8) = Dsln(8) - Dsln(7)*Dsln(2)
    Dsln(9) = Dsln(9) - Dsln(7)*Dsln(4) - Dsln(8)*Dsln(5)
    t(1) = Dsln(7)*Dsln(1)
    t(2) = Dsln(8)*Dsln(3)
    t(3) = Dsln(9)*Dsln(6)
    Dsln(10) = 1.0D0/(Dsln(10)-t(1)*Dsln(7)-t(2)*Dsln(8)-t(3)*Dsln(9))
    Dsln(7) = t(1)
    Dsln(8) = t(2)
    Dsln(9) = t(3)
    Dsln(12) = Dsln(12) - Dsln(11)*Dsln(2)
    Dsln(13) = Dsln(13) - Dsln(11)*Dsln(4) - Dsln(12)*Dsln(5)
    Dsln(14) = Dsln(14) - Dsln(11)*Dsln(7) - Dsln(12)*Dsln(8) - Dsln(13)*Dsln(9)
    t(1) = Dsln(11)*Dsln(1)
    t(2) = Dsln(12)*Dsln(3)
    t(3) = Dsln(13)*Dsln(6)
    t(4) = Dsln(14)*Dsln(10)
    Dsln(15) = 1.0D0/(Dsln(15)-t(1)*Dsln(11)-t(2)*Dsln(12)-t(3)*Dsln(13)-t(4)*Dsln(14))
    Dsln(11) = t(1)
    Dsln(12) = t(2)
    Dsln(13) = t(3)
    Dsln(14) = t(4)
    Dsln(17) = Dsln(17) - Dsln(16)*Dsln(2)
    Dsln(18) = Dsln(18) - Dsln(16)*Dsln(4) - Dsln(17)*Dsln(5)
    Dsln(19) = Dsln(19) - Dsln(16)*Dsln(7) - Dsln(17)*Dsln(8) - Dsln(18)*Dsln(9)
    Dsln(20) = Dsln(20) - Dsln(16)*Dsln(11) - Dsln(17)*Dsln(12) - Dsln(18)*Dsln(13) - Dsln(19)*Dsln(14)
    t(1) = Dsln(16)*Dsln(1)
    t(2) = Dsln(17)*Dsln(3)
    t(3) = Dsln(18)*Dsln(6)
    t(4) = Dsln(19)*Dsln(10)
    t(5) = Dsln(20)*Dsln(15)
    Dsln(21) = 1.0D0/(Dsln(21)-t(1)*Dsln(16)-t(2)*Dsln(17)-t(3) *Dsln(18)-t(4)*Dsln(19)-t(5)*Dsln(20))
    Dsln(16) = t(1)
    Dsln(17) = t(2)
    Dsln(18) = t(3)
    Dsln(19) = t(4)
    Dsln(20) = t(5)
  end subroutine INV6

  !======================================================================!
  !> @brief INV66
  !======================================================================!
  subroutine INV66(Zln,Zz,Diag)
    implicit none
    !------
    real(kind=kreal), intent(in):: Diag(21)
    real(kind=kreal), intent(in):: Zz(36)
    real(kind=kreal), intent(out):: Zln(36)
    !------
    integer(kind=kint):: i

    do i = 0, 5
      Zln(i+7) = Zz(i+7) - Zz(i+1)*Diag(2)
      Zln(i+13) = Zz(i+13) - Zz(i+1)*Diag(4) - Zln(i+7)*Diag(5)
      Zln(i+19) = Zz(i+19) - Zz(i+1)*Diag(7) - Zln(i+7)*Diag(8) - Zln(i+13)*Diag(9)
      Zln(i+25) = Zz(i+25) - Zz(i+1)*Diag(11) - Zln(i+7)*Diag(12) - Zln(i+13)*Diag(13)&
        - Zln(i+19)*Diag(14)
      Zln(i+31) = Zz(i+31) - Zz(i+1)*Diag(16) - Zln(i+7)*Diag(17) - Zln(i+13)*Diag(18)&
        - Zln(i+19)*Diag(19) - Zln(i+25)*Diag(20)
      Zln(i+1) = Zz(i+1)*Diag(1)
      Zln(i+7) = Zln(i+7)*Diag(3)
      Zln(i+13) = Zln(i+13)*Diag(6)
      Zln(i+19) = Zln(i+19)*Diag(10)
      Zln(i+25) = Zln(i+25)*Diag(15)
      Zln(i+31) = Zln(i+31)*Diag(21)
      Zln(i+25) = Zln(i+25) - Zln(i+31)*Diag(20)
      Zln(i+19) = Zln(i+19) - Zln(i+31)*Diag(19) - Zln(i+25)*Diag(14)
      Zln(i+13) = Zln(i+13) - Zln(i+31)*Diag(18) - Zln(i+25)*Diag(13)- Zln(i+19)*Diag(9)
      Zln(i+7) = Zln(i+7) - Zln(i+31)*Diag(17) - Zln(i+25)*Diag(12)- Zln(i+19)*Diag(8)&
        - Zln(i+13)*Diag(5)
      Zln(i+1) = Zln(i+1) - Zln(i+31)*Diag(16) - Zln(i+25)*Diag(11)- Zln(i+19)*Diag(7)&
        - Zln(i+13)*Diag(4) - Zln(i+7)*Diag(2)
    enddo
  end subroutine INV66

  !======================================================================!
  !> @brief INVX
  !======================================================================!
  subroutine INVX(Dsln,Ndeg,Ir)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(out):: Ir
    real(kind=kreal), intent(inout):: Dsln(*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: k
    integer(kind=kint):: k0
    integer(kind=kint):: l
    integer(kind=kint):: l0
    integer(kind=kint):: ld
    integer(kind=kint):: ll
    real(kind=kreal):: t
    real(kind=kreal):: tem

    Ir = 0
    l = 1
    Dsln(1) = 1.0D0/Dsln(1)
    do i = 2, Ndeg
      ld = 0
      l0 = l
      do j = 1, i - 1
        l = l + 1
        do k = 1, j - 1
          ld = ld + 1
          Dsln(l) = Dsln(l) - Dsln(l0+k)*Dsln(ld)
        enddo
        ld = ld + 1
      enddo
      t = 0.0D0
      k0 = 0
      ll = 0
      do k = l - i + 2, l
        ll = ll + 1
        k0 = k0 + ll
        tem = Dsln(k)*Dsln(k0)
        t = t + tem*Dsln(k)
        Dsln(k) = tem
      enddo
      l = l + 1
      Dsln(l) = Dsln(l) - t
      Dsln(l) = 1.0D0/Dsln(l)
    enddo
  end subroutine INVX

  !======================================================================!
  !> @brief INVXX
  !======================================================================!
  subroutine INVXX(Zln,Zz,Diag,Ndeg)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    real(kind=kreal), intent(in):: Diag(*)
    real(kind=kreal), intent(in):: Zz(Ndeg,Ndeg)
    real(kind=kreal), intent(out):: Zln(Ndeg,Ndeg)
    !------
    integer(kind=kint):: joc
    integer(kind=kint):: l
    integer(kind=kint):: loc1
    integer(kind=kint):: m
    integer(kind=kint):: n

    Zln = Zz
    do l = 1, Ndeg, 2
      joc = 0
      do m = 1, Ndeg - 1
        joc = joc + m
        loc1 = joc + m
        do n = m + 1, Ndeg
          Zln(l,n) = Zln(l,n) - Zln(l,m)*Diag(loc1)
          Zln(l+1,n) = Zln(l+1,n) - Zln(l+1,m)*Diag(loc1)
          loc1 = loc1 + n
        enddo
      enddo
      joc = 0
      do m = 1, Ndeg
        joc = joc + m
        Zln(l,m) = Zln(l,m)*Diag(joc)
        Zln(l+1,m) = Zln(l+1,m)*Diag(joc)
      enddo
      do n = Ndeg, 2, -1
        joc = joc - 1
        do m = n - 1, 1, -1
          Zln(l,m) = Zln(l,m) - Zln(l,n)*Diag(joc)
          Zln(l+1,m) = Zln(l+1,m) - Zln(l+1,n)*Diag(joc)
          joc = joc - 1
        enddo
      enddo
    enddo
  end subroutine INVXX

  !======================================================================!
  !> @brief DDOT performs inner product of sparse vectors
  !======================================================================!
  real(kind=kreal) function DDOT(A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(N)
    real(kind=kreal), intent(in):: B(N)
    !------
    integer(kind=kint):: i
    real(kind=kreal):: s

    s = 0.0D0
    do i = 1, N
      s = s + A(i)*B(i)
    enddo
    DDOT = s
  end function DDOT

  !======================================================================!
  !> @brief D2DOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D2DOT(T,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(4,*)
    real(kind=kreal), intent(in):: B(4,*)
    real(kind=kreal), intent(out):: T(4)
    !------
    integer(kind=kint):: jj

    T(1:4) = 0.0D0

    do jj = 1, N
      T(1) = T(1) + A(1,jj)*B(1,jj) + A(3,jj)*B(3,jj)
      T(2) = T(2) + A(2,jj)*B(1,jj) + A(4,jj)*B(3,jj)
      T(3) = T(3) + A(1,jj)*B(2,jj) + A(3,jj)*B(4,jj)
      T(4) = T(4) + A(2,jj)*B(2,jj) + A(4,jj)*B(4,jj)
    enddo
  end subroutine D2DOT

  !======================================================================!
  !> @brief D3DOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D3DOT(T,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(9,*)
    real(kind=kreal), intent(in):: B(9,*)
    real(kind=kreal), intent(out):: T(9)
    !------
    integer(kind=kint):: jj

    T(1:9) = 0.0D0
    do jj = 1, N
      T(1) = T(1) + A(1,jj)*B(1,jj) + A(4,jj)*B(4,jj) + A(7,jj)*B(7,jj)
      T(2) = T(2) + A(2,jj)*B(1,jj) + A(5,jj)*B(4,jj) + A(8,jj)*B(7,jj)
      T(3) = T(3) + A(3,jj)*B(1,jj) + A(6,jj)*B(4,jj) + A(9,jj)*B(7,jj)

      T(4) = T(4) + A(1,jj)*B(2,jj) + A(4,jj)*B(5,jj) + A(7,jj)*B(8,jj)
      T(5) = T(5) + A(2,jj)*B(2,jj) + A(5,jj)*B(5,jj) + A(8,jj)*B(8,jj)
      T(6) = T(6) + A(3,jj)*B(2,jj) + A(6,jj)*B(5,jj) + A(9,jj)*B(8,jj)

      T(7) = T(7) + A(1,jj)*B(3,jj) + A(4,jj)*B(6,jj) + A(7,jj)*B(9,jj)
      T(8) = T(8) + A(2,jj)*B(3,jj) + A(5,jj)*B(6,jj) + A(8,jj)*B(9,jj)
      T(9) = T(9) + A(3,jj)*B(3,jj) + A(6,jj)*B(6,jj) + A(9,jj)*B(9,jj)
    enddo
  end subroutine D3DOT

  !======================================================================!
  !> @brief D3DOTL performs inner product of sparse vectors
  !======================================================================!
  subroutine D3DOTL(T,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(9,*)
    real(kind=kreal), intent(in):: B(9,*)
    real(kind=kreal), intent(out):: T(6)
    !------
    integer(kind=kint):: jj

    T(1:6) = 0.0D0
    do jj = 1, N
      T(1) = T(1) + A(1,jj)*B(1,jj) + A(4,jj)*B(4,jj) + A(7,jj)*B(7,jj)
      T(2) = T(2) + A(2,jj)*B(1,jj) + A(5,jj)*B(4,jj) + A(8,jj)*B(7,jj)

      T(3) = T(3) + A(2,jj)*B(2,jj) + A(5,jj)*B(5,jj) + A(8,jj)*B(8,jj)
      T(4) = T(4) + A(3,jj)*B(1,jj) + A(6,jj)*B(4,jj) + A(9,jj)*B(7,jj)

      T(5) = T(5) + A(3,jj)*B(2,jj) + A(6,jj)*B(5,jj) + A(9,jj)*B(8,jj)
      T(6) = T(6) + A(3,jj)*B(3,jj) + A(6,jj)*B(6,jj) + A(9,jj)*B(9,jj)
    enddo
  end subroutine D3DOTL

  !======================================================================!
  !> @brief DXDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine DXDOT(Ndeg,T,A,B,L)
    implicit none
    !------
    integer(kind=kint), intent(in):: L
    integer(kind=kint), intent(in):: Ndeg
    real(kind=kreal), intent(in):: A(Ndeg,Ndeg,*)
    real(kind=kreal), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=kreal), intent(out):: T(Ndeg,Ndeg)
    !------
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: m
    integer(kind=kint):: n

    do n = 1, Ndeg
      do m = 1, Ndeg
        T(n,m) = 0.0D0
        do k = 1, Ndeg
          do jj = 1, L
            T(n,m) = T(n,m) + A(n,k,jj)*B(m,k,jj)
          enddo
        enddo
      enddo
    enddo
  end subroutine DXDOT

  !======================================================================!
  !> @brief DXDOTL performs inner product of sparse vectors
  !======================================================================!
  subroutine DXDOTL(Ndeg,T,A,B,L)
    implicit none
    !------
    integer(kind=kint), intent(in):: L
    integer(kind=kint), intent(in):: Ndeg
    real(kind=kreal), intent(in):: A(Ndeg,Ndeg,*)
    real(kind=kreal), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=kreal), intent(out):: T(Ndeg,Ndeg)
    !------
    integer(kind=kint):: jj
    integer(kind=kint):: k
    integer(kind=kint):: m
    integer(kind=kint):: n

    do n = 1, Ndeg
      do m = 1, n
        T(n,m) = 0.0D0
        do k = 1, Ndeg
          do jj = 1, L
            T(n,m) = T(n,m) + A(n,k,jj)*B(m,k,jj)
          enddo
        enddo
      enddo
    enddo
  end subroutine DXDOTL

  !======================================================================!
  !> @brief VPROD
  !======================================================================!
  subroutine VPROD(A,B,C,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(N)
    real(kind=kreal), intent(in):: B(N)
    real(kind=kreal), intent(out):: C(N)
    !------
    C(1:N) = A(1:N)*B(1:N)
  end subroutine VPROD

  !======================================================================!
  !> @brief V2PROD
  !======================================================================!
  subroutine V2PROD(A,B,C,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(4,N)
    real(kind=kreal), intent(in):: B(3,N)
    real(kind=kreal), intent(out):: C(4,N)
    !------
    integer(kind=kint):: i

    do i = 1, N
      C(3,i) = A(3,i) - A(1,i)*B(2,i)
      C(1,i) = A(1,i)*B(1,i)
      C(3,i) = C(3,i)*B(3,i)
      C(1,i) = C(1,i) - C(3,i)*B(2,i)

      C(4,i) = A(4,i) - A(2,i)*B(2,i)
      C(2,i) = A(2,i)*B(1,i)
      C(4,i) = C(4,i)*B(3,i)
      C(2,i) = C(2,i) - C(4,i)*B(2,i)
    enddo
  end subroutine V2PROD

  !======================================================================!
  !> @brief V3PROD
  !======================================================================!
  subroutine V3PROD(Zln,Diag,Zz,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: Diag(6,N)
    real(kind=kreal), intent(in):: Zln(9,N)
    real(kind=kreal), intent(out):: Zz(9,N)
    !------
    integer(kind=kint):: i

    do i = 1, N
      Zz(4,i) = Zln(4,i) - Zln(1,i)*Diag(2,i)
      Zz(7,i) = Zln(7,i) - Zln(1,i)*Diag(4,i) - Zz(4,i)*Diag(5,i)
      Zz(1,i) = Zln(1,i)*Diag(1,i)
      Zz(4,i) = Zz(4,i)*Diag(3,i)
      Zz(7,i) = Zz(7,i)*Diag(6,i)
      Zz(4,i) = Zz(4,i) - Zz(7,i)*Diag(5,i)
      Zz(1,i) = Zz(1,i) - Zz(4,i)*Diag(2,i) - Zz(7,i)*Diag(4,i)

      Zz(5,i) = Zln(5,i) - Zln(2,i)*Diag(2,i)
      Zz(8,i) = Zln(8,i) - Zln(2,i)*Diag(4,i) - Zz(5,i)*Diag(5,i)
      Zz(2,i) = Zln(2,i)*Diag(1,i)
      Zz(5,i) = Zz(5,i)*Diag(3,i)
      Zz(8,i) = Zz(8,i)*Diag(6,i)
      Zz(5,i) = Zz(5,i) - Zz(8,i)*Diag(5,i)
      Zz(2,i) = Zz(2,i) - Zz(5,i)*Diag(2,i) - Zz(8,i)*Diag(4,i)

      Zz(6,i) = Zln(6,i) - Zln(3,i)*Diag(2,i)
      Zz(9,i) = Zln(9,i) - Zln(3,i)*Diag(4,i) - Zz(6,i)*Diag(5,i)
      Zz(3,i) = Zln(3,i)*Diag(1,i)
      Zz(6,i) = Zz(6,i)*Diag(3,i)
      Zz(9,i) = Zz(9,i)*Diag(6,i)
      Zz(6,i) = Zz(6,i) - Zz(9,i)*Diag(5,i)
      Zz(3,i) = Zz(3,i) - Zz(6,i)*Diag(2,i) - Zz(9,i)*Diag(4,i)
    enddo
  end subroutine V3PROD

  !======================================================================!
  !> @brief VXPROD
  !======================================================================!
  subroutine VXPROD(Ndeg,Ndegl,Zln,Diag,Zz,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    real(kind=kreal), intent(in):: Diag(Ndegl,N)
    real(kind=kreal), intent(in):: Zln(Ndeg*Ndeg,N)
    real(kind=kreal), intent(out):: Zz(Ndeg*Ndeg,N)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: N

    do i = 1, N
      call INVXX(Zz(1,i),Zln(1,i),Diag(1,i),Ndeg)
    enddo
  end subroutine VXPROD

  !======================================================================!
  !> @brief VCOPY
  !======================================================================!
  subroutine VCOPY(A,C,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(N)
    real(kind=kreal), intent(out):: C(N)
    !------
    C = A
  end subroutine VCOPY

  !======================================================================!
  !> @brief NUSOL0 performs forward elimination and backward substitution
  !     (i/o)
  !           r_h_s    on entry     right hand side vector
  !                    on exit      solution vector
  !           iv       communication array
  !  updates STAge
  !======================================================================!
  subroutine NUSOL0(R_h_s,FCT,Ir)
    implicit none
    !------
    real(kind=kreal), intent(inout):: R_h_s(:)
    type (cholesky_factor), intent(inout):: FCT
    integer(kind=kint), intent(out):: Ir
    !------
    integer(kind=kint):: ndegl
    integer(kind=kint):: ierror
    real(kind=kreal), pointer :: wk(:)

    if ( FCT%STAge/=30 .and. FCT%STAge/=40 ) then
      Ir = 50
      return
    else
      Ir = 0
    endif

    allocate (wk(FCT%NDEg*FCT%NEQns),stat=IERror)
    if ( IERror/=0 ) then
      write (*,*) "##Error: not enough memory"
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
    endif
    !rmiv
    ndegl = FCT%NDEg*(FCT%NDEg+1)
    ndegl = ndegl/2

    select case( FCT%NDEg )
      case (1)
        call NUSOL1(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,FCT%IPErm,R_h_s,wk,FCT%NEQns,FCT%NSTop)
      case (2)
        call NUSOL2(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,FCT%IPErm,R_h_s,wk,FCT%NEQns,FCT%NSTop)
      case (3)
        call NUSOL3(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,FCT%IPErm,R_h_s,wk,FCT%NEQns,FCT%NSTop)
      case (6)
        call NUSOLX(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,FCT%IPErm,R_h_s,wk,FCT%NEQns,FCT%NSTop,FCT%NDEg,ndegl)
      case default
        call NUSOLX(FCT%XLNzr,FCT%COLno,FCT%DSLn,FCT%ZLN,FCT%DIAg,FCT%IPErm,R_h_s,wk,FCT%NEQns,FCT%NSTop,FCT%NDEg,ndegl)
    endselect

    FCT%STAge = 40
    deallocate (wk)
  end subroutine NUSOL0

  !======================================================================!
  !> @brief NUSOL1 performs forward elimination and backward substitution
  !======================================================================!
  subroutine NUSOL1(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Iperm(:)
    real(kind=kreal), intent(in):: Zln(:)
    real(kind=kreal), intent(in):: Diag(:)
    real(kind=kreal), intent(in):: Dsln(:)
    real(kind=kreal), intent(inout):: B(:)
    real(kind=kreal), intent(out):: Wk(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    ! forward
    do i = 1, Neqns
      Wk(i) = B(Iperm(i))
    enddo
    joc = 1
    do i = 1, Neqns
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) Wk(i) = Wk(i) - SPDOT2(Wk,Zln,Colno,ks,ke)
      if ( i>Nstop ) then
        Wk(i) = Wk(i) - DDOT(Wk(Nstop:),Dsln(joc:),i-Nstop)
        joc = joc + i - Nstop
      endif
    enddo
    do i = 1, Neqns
      Wk(i) = Wk(i)*Diag(i)
    enddo
    ! back ward
    do i = Neqns, 1, -1
      if ( i>=Nstop ) then
        do j = i - 1, Nstop, -1
          joc = joc - 1
          Wk(j) = Wk(j) - Wk(i)*Dsln(joc)
        enddo
      endif
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) then
        do k = ks, ke
          j = Colno(k)
          Wk(j) = Wk(j) - Wk(i)*Zln(k)
        enddo
      endif
    enddo
    ! permutation
    do i = 1, Neqns
      B(Iperm(i)) = Wk(i)
    enddo
  end subroutine NUSOL1

  !======================================================================!
  !> @brief NUSOL2 performs forward elimination and backward substitution
  !======================================================================!
  subroutine NUSOL2(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Iperm(:)
    real(kind=kreal), intent(in):: Zln(4,*)
    real(kind=kreal), intent(in):: Diag(3,*)
    real(kind=kreal), intent(in):: Dsln(4,*)
    real(kind=kreal), intent(inout):: B(2,*)
    real(kind=kreal), intent(out):: Wk(2,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    ! forward
    do i = 1, Neqns
      Wk(1,i) = B(1,Iperm(i))
      Wk(2,i) = B(2,Iperm(i))
    enddo
    joc = 1
    do i = 1, Neqns
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) call S2PDOT(Wk(1,i),Wk,Zln,Colno,ks,ke)
      if ( i>Nstop ) then
        call D2SDOT(Wk(1,i),Wk(1,Nstop),Dsln(1,joc),i-Nstop)
        joc = joc + i - Nstop
      endif
    enddo
    do i = 1, Neqns
      Wk(2,i) = Wk(2,i) - Wk(1,i)*Diag(2,i)
      Wk(1,i) = Wk(1,i)*Diag(1,i)
      Wk(2,i) = Wk(2,i)*Diag(3,i)
      Wk(1,i) = Wk(1,i) - Wk(2,i)*Diag(2,i)
    enddo
    ! back ward
    do i = Neqns, 1, -1
      if ( i>=Nstop ) then
        do j = i - 1, Nstop, -1
          joc = joc - 1
          Wk(1,j) = Wk(1,j) - Wk(1,i)*Dsln(1,joc) - Wk(2,i)*Dsln(2,joc)
          Wk(2,j) = Wk(2,j) - Wk(1,i)*Dsln(3,joc) - Wk(2,i)*Dsln(4,joc)
        enddo
      endif
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) then
        do k = ks, ke
          j = Colno(k)
          Wk(1,j) = Wk(1,j) - Wk(1,i)*Zln(1,k) - Wk(2,i)*Zln(2,k)
          Wk(2,j) = Wk(2,j) - Wk(1,i)*Zln(3,k) - Wk(2,i)*Zln(4,k)
        enddo
      endif
    enddo
    ! permutation
    do i = 1, Neqns
      B(1,Iperm(i)) = Wk(1,i)
      B(2,Iperm(i)) = Wk(2,i)
    enddo
  end subroutine NUSOL2

  !======================================================================!
  !> @brief NUSOL3 performs forward elimination and backward substitution
  !======================================================================!
  subroutine NUSOL3(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Iperm(:)
    real(kind=kreal), intent(in):: Zln(9,*)
    real(kind=kreal), intent(in):: Diag(6,*)
    real(kind=kreal), intent(in):: Dsln(9,*)
    real(kind=kreal), intent(inout):: B(3,*)
    real(kind=kreal), intent(out):: Wk(3,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks

    ! forward
    do i = 1, Neqns
      Wk(1,i) = B(1,Iperm(i))
      Wk(2,i) = B(2,Iperm(i))
      Wk(3,i) = B(3,Iperm(i))
    enddo
    joc = 1
    do i = 1, Neqns
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) call S3PDOT(Wk(1,i),Wk,Zln,Colno,ks,ke)
      if ( i>Nstop ) then
        call D3SDOT(Wk(1,i),Wk(1,Nstop),Dsln(1,joc),i-Nstop)
        joc = joc + i - Nstop
      endif
    enddo
    do i = 1, Neqns
      Wk(2,i) = Wk(2,i) - Wk(1,i)*Diag(2,i)
      Wk(3,i) = Wk(3,i) - Wk(1,i)*Diag(4,i) - Wk(2,i)*Diag(5,i)
      Wk(1,i) = Wk(1,i)*Diag(1,i)
      Wk(2,i) = Wk(2,i)*Diag(3,i)
      Wk(3,i) = Wk(3,i)*Diag(6,i)
      Wk(2,i) = Wk(2,i) - Wk(3,i)*Diag(5,i)
      Wk(1,i) = Wk(1,i) - Wk(2,i)*Diag(2,i) - Wk(3,i)*Diag(4,i)
    enddo
    ! back ward
    do i = Neqns, 1, -1
      if ( i>=Nstop ) then
        do j = i - 1, Nstop, -1
          joc = joc - 1
          Wk(1,j) = Wk(1,j) - Wk(1,i)*Dsln(1,joc) - Wk(2,i)*Dsln(2,joc) - Wk(3,i)*Dsln(3,joc)
          Wk(2,j) = Wk(2,j) - Wk(1,i)*Dsln(4,joc) - Wk(2,i)*Dsln(5,joc) - Wk(3,i)*Dsln(6,joc)
          Wk(3,j) = Wk(3,j) - Wk(1,i)*Dsln(7,joc) - Wk(2,i)*Dsln(8,joc) - Wk(3,i)*Dsln(9,joc)
        enddo
      endif
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) then
        do k = ks, ke
          j = Colno(k)
          Wk(1,j) = Wk(1,j) - Wk(1,i)*Zln(1,k) - Wk(2,i)*Zln(2,k) - Wk(3,i)*Zln(3,k)
          Wk(2,j) = Wk(2,j) - Wk(1,i)*Zln(4,k) - Wk(2,i)*Zln(5,k) - Wk(3,i)*Zln(6,k)
          Wk(3,j) = Wk(3,j) - Wk(1,i)*Zln(7,k) - Wk(2,i)*Zln(8,k) - Wk(3,i)*Zln(9,k)
        enddo
      endif
    enddo
    ! permutation
    do i = 1, Neqns
      B(1,Iperm(i)) = Wk(1,i)
      B(2,Iperm(i)) = Wk(2,i)
      B(3,Iperm(i)) = Wk(3,i)
    enddo
  end subroutine NUSOL3

  !======================================================================!
  !> @brief NUSOLX performs forward elimination and backward substitution
  !======================================================================!
  subroutine NUSOLX(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop,Ndeg,Ndegl)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Ndegl
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nstop
    integer(kind=kint), intent(in):: Xlnzr(:)
    integer(kind=kint), intent(in):: Colno(:)
    integer(kind=kint), intent(in):: Iperm(:)
    real(kind=kreal), intent(in):: Zln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(in):: Diag(Ndegl,*)
    real(kind=kreal), intent(in):: Dsln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: B(Ndeg,*)
    real(kind=kreal), intent(out):: Wk(Ndeg,*)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: j
    integer(kind=kint):: joc
    integer(kind=kint):: joc1
    integer(kind=kint):: k
    integer(kind=kint):: ke
    integer(kind=kint):: ks
    integer(kind=kint):: l
    integer(kind=kint):: loc1
    integer(kind=kint):: locd
    integer(kind=kint):: m
    integer(kind=kint):: n

    ! forward
    do l = 1, Ndeg
      Wk(l,1:Neqns) = B(l,Iperm(1:Neqns))
    enddo
    joc = 1
    do i = 1, Neqns
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) call SXPDOT(Ndeg,Wk(1,i),Wk,Zln,Colno,ks,ke)
      if ( i>Nstop ) then
        joc1 = i - Nstop
        call DXSDOT(Ndeg,Wk(1,i),Wk(1,Nstop),Dsln(1,1,joc),joc1)
        joc = joc + joc1
      endif
    enddo
    do i = 1, Neqns
      locd = 0
      do m = 1, Ndeg - 1
        locd = locd + m
        loc1 = locd + m
        do n = m + 1, Ndeg
          Wk(n,i) = Wk(n,i) - Wk(m,i)*Diag(loc1,i)
          loc1 = loc1 + n
        enddo
      enddo
      locd = 0
      do m = 1, Ndeg
        locd = locd + m
        Wk(m,i) = Wk(m,i)*Diag(locd,i)
      enddo
      do n = Ndeg, 2, -1
        locd = locd - 1
        do m = n - 1, 1, -1
          Wk(m,i) = Wk(m,i) - Wk(n,i)*Diag(locd,i)
          locd = locd - 1
        enddo
      enddo
    enddo
    ! back ward
    do i = Neqns, 1, -1
      if ( i>=Nstop ) then
        do j = i - 1, Nstop, -1
          joc = joc - 1
          do m = 1, Ndeg
            do n = 1, Ndeg
              Wk(m,j) = Wk(m,j) - Wk(n,i)*Dsln(n,m,joc)
            enddo
          enddo
        enddo
      endif
      ks = Xlnzr(i)
      ke = Xlnzr(i+1) - 1
      if ( ke>=ks ) then
        do k = ks, ke
          j = Colno(k)
          do m = 1, Ndeg
            do n = 1, Ndeg
              Wk(m,j) = Wk(m,j) - Wk(n,i)*Zln(n,m,k)
            enddo
          enddo
        enddo
      endif
    enddo
    ! permutation
    do l = 1, Ndeg
      B(l,Iperm(1:Neqns)) = Wk(l,1:Neqns)
    enddo
  end subroutine NUSOLX

  !======================================================================!
  !> @brief SPDOT2 performs inner product of sparse vectors
  !======================================================================!
  real(kind=kreal) function SPDOT2(B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ke
    integer(kind=kint), intent(in):: Ks
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Zln(:)
    real(kind=kreal), intent(in):: B(:)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jj
    real(kind=kreal):: s

    s = 0.0D0
    do jj = Ks, Ke
      j = Colno(jj)
      s = s + Zln(jj)*B(j)
    enddo
    SPDOT2 = s
  end function SPDOT2

  !======================================================================!
  !> @brief S2PDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine S2PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ke
    integer(kind=kint), intent(in):: Ks
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Zln(4,*)
    real(kind=kreal), intent(in):: B(2,*)
    real(kind=kreal), intent(inout):: Bi(2)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jj

    do jj = Ks, Ke
      j = Colno(jj)
      Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(3,jj)*B(2,j)
      Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(4,jj)*B(2,j)
    enddo
  end subroutine S2PDOT

  !======================================================================!
  !> @brief S3PDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine S3PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ke
    integer(kind=kint), intent(in):: Ks
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Zln(9,*)
    real(kind=kreal), intent(in):: B(3,*)
    real(kind=kreal), intent(inout):: Bi(3)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jj

    do jj = Ks, Ke
      j = Colno(jj)
      Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(4,jj)*B(2,j) - Zln(7,jj)*B(3,j)
      Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(5,jj)*B(2,j) - Zln(8,jj)*B(3,j)
      Bi(3) = Bi(3) - Zln(3,jj)*B(1,j) - Zln(6,jj)*B(2,j) - Zln(9,jj)*B(3,j)
    enddo
  end subroutine S3PDOT

  !======================================================================!
  !> @brief S6PDOT  performs inner product of sparse vectors
  !======================================================================!
  subroutine S6PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ke
    integer(kind=kint), intent(in):: Ks
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Zln(36,*)
    real(kind=kreal), intent(in):: B(6,*)
    real(kind=kreal), intent(inout):: Bi(6)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jj

    do jj = Ks, Ke
      j = Colno(jj)
      Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(7,jj)*B(2,j) - Zln(13,jj)*B(3,j)&
        - Zln(19,jj)*B(4,j) - Zln(25,jj)*B(5,j) - Zln(31,jj)*B(6,j)
      Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(8,jj)*B(2,j) - Zln(14,jj)*B(3,j)&
        - Zln(20,jj)*B(4,j) - Zln(26,jj)*B(5,j) - Zln(32,jj)*B(6,j)
      Bi(3) = Bi(3) - Zln(3,jj)*B(1,j) - Zln(9,jj)*B(2,j) - Zln(15,jj)*B(3,j)&
        - Zln(21,jj)*B(4,j) - Zln(27,jj)*B(5,j) - Zln(33,jj)*B(6,j)
      Bi(4) = Bi(4) - Zln(4,jj)*B(1,j) - Zln(10,jj)*B(2,j) - Zln(16,jj)*B(3,j)&
        - Zln(22,jj)*B(4,j) - Zln(28,jj)*B(5,j) - Zln(34,jj)*B(6,j)
      Bi(5) = Bi(5) - Zln(5,jj)*B(1,j) - Zln(11,jj)*B(2,j) - Zln(17,jj)*B(3,j)&
        - Zln(23,jj)*B(4,j) - Zln(29,jj)*B(5,j) - Zln(35,jj)*B(6,j)
      Bi(6) = Bi(6) - Zln(6,jj)*B(1,j) - Zln(12,jj)*B(2,j) - Zln(18,jj)*B(3,j)&
        - Zln(25,jj)*B(4,j) - Zln(30,jj)*B(5,j) - Zln(36,jj)*B(6,j)
    enddo
  end subroutine S6PDOT

  !======================================================================!
  !> @brief SXPDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine SXPDOT(Ndeg,Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ke
    integer(kind=kint), intent(in):: Ks
    integer(kind=kint), intent(in):: Ndeg
    integer(kind=kint), intent(in):: Colno(:)
    real(kind=kreal), intent(in):: Zln(Ndeg,Ndeg,*)
    real(kind=kreal), intent(in):: B(Ndeg,*)
    real(kind=kreal), intent(inout):: Bi(Ndeg)
    !------
    integer(kind=kint):: j
    integer(kind=kint):: jj
    integer(kind=kint):: m
    integer(kind=kint):: n

    do jj = Ks, Ke
      j = Colno(jj)
      do m = 1, Ndeg
        do n = 1, Ndeg
          Bi(n) = Bi(n) - Zln(n,m,jj)*B(m,j)
        enddo
      enddo
    enddo
  end subroutine SXPDOT

  !======================================================================!
  !> @brief D2SDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D2SDOT(Wi,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(2,*)
    real(kind=kreal), intent(in):: B(4,*)
    real(kind=kreal), intent(inout):: Wi(2)
    !------
    integer(kind=kint):: jj

    do jj = 1, N
      Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(3,jj)
      Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(4,jj)
    enddo
  end subroutine D2SDOT

  !======================================================================!
  !> @brief D3SDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D3SDOT(Wi,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: N
    real(kind=kreal), intent(in):: A(3,*)
    real(kind=kreal), intent(in):: B(9,*)
    real(kind=kreal), intent(inout):: Wi(3)
    !------
    integer(kind=kint):: jj

    do jj = 1, N
      Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(4,jj) - A(3,jj)*B(7,jj)
      Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(5,jj) - A(3,jj)*B(8,jj)
      Wi(3) = Wi(3) - A(1,jj)*B(3,jj) - A(2,jj)*B(6,jj) - A(3,jj)*B(9,jj)
    enddo
  end subroutine D3SDOT

  !======================================================================!
  !> @brief DXSDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine DXSDOT(Ndeg,Wi,A,B,N)
    implicit none
    !------
    integer(kind=kint), intent(in):: Ndeg
    real(kind=kreal), intent(in):: A(Ndeg,*)
    real(kind=kreal), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=kreal), intent(inout):: Wi(Ndeg)
    integer(kind=kint), intent(inout):: N
    !------
    integer(kind=kint):: jj
    integer(kind=kint):: m

    do jj = 1, N
      do m = 1, Ndeg
        do N = 1, Ndeg
          Wi(N) = Wi(N) - B(N,m,jj)*A(m,jj)
        enddo
      enddo
    enddo
  end subroutine DXSDOT

end module HECMW_SOLVER_DIRECT
