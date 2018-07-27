!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_SOLVE_DIRECT is a program for the matrix direct solver
!         SOLVER=DIRECT
!
!> @author coded by t.arakawa of RIST on 040316
!>         Modified G. Prabhakar, RIST, June 7 2004
!----------------------------------------------------------------------
module HECMW_SOLVER_DIRECT
  implicit none

  integer(kind=4), private :: len_colno
  integer(kind=4), private :: nstop
  integer(kind=4), private :: stage
  integer(kind=4), private :: neqns
  integer(kind=4), private :: nttbr
  integer(kind=4), private :: isym
  integer(kind=4), private :: ndeg
  integer(kind=4), private :: irr
  integer(kind=4), private :: len_dsln
  integer(kind=4), private :: len_iv
  integer(kind=4), private :: total

  integer(kind=4), private, dimension(:), pointer :: jcol
  integer(kind=4), private, dimension(:), pointer :: irow
  integer(kind=4), private, dimension(:), pointer :: zpiv
  integer(kind=4), private, dimension(:), pointer :: iperm
  integer(kind=4), private, dimension(:), pointer :: invp
  integer(kind=4), private, dimension(:), pointer :: parent
  integer(kind=4), private, dimension(:), pointer :: nch
  integer(kind=4), private, dimension(:), pointer :: xlnzr
  integer(kind=4), private, dimension(:), pointer :: colno
  !*Work arrays
  integer(kind=4), private, dimension(:), pointer :: jcpt
  integer(kind=4), private, dimension(:), pointer :: jcolno
  integer(kind=4), private, dimension(:), pointer :: ia
  integer(kind=4), private, dimension(:), pointer :: ja
  integer(kind=4), private, dimension(:), pointer :: deg
  integer(kind=4), private, dimension(:), pointer :: marker
  integer(kind=4), private, dimension(:), pointer :: rchset
  integer(kind=4), private, dimension(:), pointer :: nbrhd
  integer(kind=4), private, dimension(:), pointer :: qsize
  integer(kind=4), private, dimension(:), pointer :: qlink
  integer(kind=4), private :: nofsub
  integer(kind=4), private, dimension(:), pointer :: adjncy
  integer(kind=4), private, dimension(:), pointer :: btree
  integer(kind=4), private, dimension(:), pointer :: pordr
  integer(kind=4), private, dimension(:), pointer :: adjncp
  integer(kind=4), private, dimension(:), pointer :: xleaf
  integer(kind=4), private, dimension(:), pointer :: leaf
  integer(kind=4), private, dimension(:), pointer :: indx

  real(kind=8), private, dimension(:), pointer :: val
  real(kind=8), private, dimension(:), pointer :: temp
  real(kind=8), private, dimension(:), pointer :: diag
  real(kind=8), private, dimension(:), pointer :: zln
  real(kind=8), private, dimension(:), pointer :: dsln
  !*Timing
  real(kind=8), private, dimension(10) :: tom
  !*Allocation variables
  integer(kind=4), private :: ialoc
  integer(kind=4), private :: raloc
  integer(kind=4), private :: ierror
  !*MCHDPN
  integer:: LRAtio
  real(kind=8), private:: EPSm
  real(kind=8), private:: RMAx
  real(kind=8), private:: RMIn
  !*QAZ
  integer, private:: ISEed
  integer, private:: ISEm
  integer, private:: IXXx
  !*DEBUG
  integer, private:: IDBg
  integer, private:: IDBg1

contains
  !----------------------------------------------------------------------
  !> @brief HECMW_SOLVE_DIRECT is a program for the matrix solver
  !     and is also available for performance tests
  !     to solve Ax=b for x, nozero pattern and values must be give.
  !   arrays
  !     jcol     column entrys
  !     irow     row entrys
  !     val      its value
  !     v        interface array used through matini, staijx, nufctx,
  !              nusolx
  !----------------------------------------------------------------------
  subroutine HECMW_SOLVE_DIRECT(hecMESH,hecMAT,Ifmsg)
    use HECMW_UTIL
    use HECMW_MATRIX_ASS
    use HECMW_MATRIX_DUMP
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    integer, intent(in):: Ifmsg
    type (HECMWST_MATRIX), intent(out)::hecMAT
    !------
    integer:: i98
    integer:: i97
    integer:: ir
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5

    RMAx = 8.988D+307
    RMIn = 4.941D-300
    EPSm = 2.220D-16
    LRAtio = 2
    ISEed = 1
    ir = 0

    call HECMW_MAT_DUMP(hecMAT,hecMESH)

    call PTIME(t1)

    !*EHM HECMW June 7 2004
    i98 = hecMAT%IARRAY(98)
    if ( hecMAT%IARRAY(98)==1 ) then
      !* Interface to symbolic factorization
      call SETIJ(hecMESH,hecMAT)

      !* Symbolic factorization
      call MATINI(ir)
      hecMAT%IARRAY(98) = 0
      write (6,*) "symbolic fct done"
    endif
    call PTIME(t2)
    t3 = t2

    i97 = hecMAT%IARRAY(97)
    if ( hecMAT%IARRAY(97)==1 ) then
      !* Interface to numeric factorization
      call NUFORM(hecMESH,hecMAT,ir)
      call PTIME(t3)

      !* Numeric factorization
      call NUFCT0(ir)
      hecMAT%IARRAY(97) = 0

      !*Memory Details
      write (*,*) '*-----------------------------------*'
      write (*,*) '|   Direct  Solver  Memory  Usage   |'
      write (*,*) '*-----------------------------------*'
      write (*,*) 'INTEGER memory: ', real(IALoc*4)/real(1048576), 'MB'
      write (*,*) 'REAL*8  memory: ', real(RALoc*8)/real(1048576), 'MB'
      write (*,*) 'TOTAL   memory: ', real((RALoc*2+IALoc)*4)/real(1048576), 'MB'
      write (*,*) '*-----------------------------------*'
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
      !WINDEBUG
      write (Ifmsg,*) 'ERROR in nufct0. ir = ', ir
      stop
    endif

    TOM(1) = t2-t1
    TOM(2) = t3-t2
    TOM(3) = t4-t3

    !* Solve
    !* Backsubstitute
    call NUSOL0(hecMAT%B,ir)
    call PTIME(t5)
    !* Errors 4
    if ( ir/=0 ) then
      !WINDEBUG
      write (Ifmsg,*) 'error in nusol0. irr = ', ir
      stop
    endif
    call HECMW_MAT_DUMP_SOLUTION(hecMAT)
  end subroutine HECMW_SOLVE_DIRECT

  !======================================================================!
  !> @brief PTIME
  !======================================================================!
  subroutine PTIME(Cputim)
    use HECMW_UTIL
    implicit none
    !------
    real(kind=8), intent(out):: Cputim
    !------
    ! cpu time by hour
    Cputim = HECMW_WTIME()
  end subroutine PTIME

  !======================================================================!
  !> @brief SETIJ
  !======================================================================!
  subroutine SETIJ(hecMESH,hecMAT)
    use HECMW_UTIL
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    type (HECMWST_MATRIX), intent(in)::hecMAT
    !------
    integer:: i
    integer:: ierr
    integer:: j
    integer:: k
    integer:: kk
    integer:: ntotal
    integer(kind=kint):: numnp
    integer(kind=kint):: ndof
    integer(kind=kint):: ndof2

    numnp = hecMAT%NP
    ndof = hecMESH%N_DOF
    ntotal = numnp*ndof

    !*NUFACT variables
    NEQns = numnp
    NDEg = ndof
    NTTbr = hecMAT%NP + hecMAT%NPL
    !+hecMAT%NPU if unsymmetric
    ISYm = 0

    !*Allocations
    allocate (IROw(NTTbr),stat=ierr)
    allocate (JCOl(NTTbr),stat=ierr)
    if ( ierr/=0 ) stop "Allocation error: irow/jcol"

    kk = 0
    ndof2 = ndof*ndof
    do j = 1, numnp
      !*Diagonal
      kk = kk + 1
      IROw(kk) = j
      JCOl(kk) = j
      !*Lower
      do k = hecMAT%INDEXL(j-1) + 1, hecMAT%INDEXL(j)
        i = hecMAT%ITEML(k)
        kk = kk + 1
        IROw(kk) = j
        JCOl(kk) = i
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
  !        iv        comunication array. v is the original name
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
  !======================================================================!
  subroutine MATINI(Ir)
    implicit none
    !------
    integer, intent(out):: Ir
    !------
    integer:: ir1
    integer:: iv1
    integer:: izz
    integer:: izz0
    integer:: ladp
    integer:: last
    integer:: lbtree
    integer:: lcolno
    integer:: lcpt
    integer:: left
    integer:: lenv
    integer:: lenv2
    integer:: lia
    integer:: lja
    integer:: lleaf
    integer:: lncol
    integer:: lnleaf
    integer:: lpordr
    integer:: lwk1
    integer:: lwk2
    integer:: lwk3
    integer:: lwk4
    integer:: lwk5
    integer:: lwk6
    integer:: lwk7
    integer:: lwk8
    integer:: lxleaf
    integer:: maxl
    integer:: neqns1
    integer:: neqnsz

    IDBg = 0
    izz0 = 0

    Ir = 0
    lenv = 1000000000
    lenv2 = LRAtio*lenv
    neqns1 = NEQns + 2
    LEN_dsln = lenv2
    LEN_iv = lenv2
    iv1 = 51

    !*Initialize allocation measure variables
    IALoc = 0
    RALoc = 0
    !
    !  set z pivot
    !
    allocate (ZPIv(NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, zpiv: SUB. matini"
    call ZPIVOT(NEQns,neqnsz,NTTbr,JCOl,IROw,ZPIv,ir1)
    if ( ir1/=0 ) then
      Ir = ir1
      return
    endif

    !
    !  build jcpt,jcolno
    !
    !rmiv
    lcpt = iv1 + neqns1
    lcolno = lcpt + 2*NTTbr
    left = lcolno + 2*NTTbr
    last = lenv2
    !rmem
    allocate (JCPt(2*NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, jcpt: SUB. matini"
    allocate (JCOlno(2*NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, jcolno: SUB. matini"
    call STSMAT(NEQns,NTTbr,IROw,JCOl,JCPt,JCOlno)
    !
    !  build ia,ja
    !
    lia = last - neqns1
    lja = lia - NTTbr*2
    last = lja
    !rmem
    allocate (IA(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, ia: SUB. matini"
    !WINDEBUG
    allocate (JA(2*NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, ja: SUB. matini"
    call STIAJA(NEQns,IA,JA,JCPt,JCOlno)

    !*Deallocation of work array
    deallocate (JCPt)
    deallocate (JCOlno)
    !
    !  get permutation vector iperm,invp
    !
    lwk1 = lja - neqns1
    lwk2 = lwk1 - neqns1
    lwk3 = lwk2 - neqns1
    lwk4 = lwk3 - neqns1
    lwk5 = lwk4 - neqns1
    lwk6 = lwk5 - neqns1
    lwk7 = lwk6 - neqns1
    lwk8 = lwk7 - 2*NTTbr
    last = lwk8

    left = iv1 + 5*neqns1

    allocate (IPErm(NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, iperm: SUB. matini"
    allocate (INVp(NEQns),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, invp: SUB. matini"
    allocate (DEG(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, deg: SUB. matini"
    allocate (MARker(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, marker: SUB. matini"
    allocate (RCHset(0:NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, rchset: SUB. matini"
    allocate (NBRhd(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, nbrhd: SUB. matini"
    allocate (QSIze(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, qsize: SUB. matini"
    allocate (QLInk(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, qlink: SUB. matini"
    allocate (ADJncy(2*NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, adjncy: SUB. matini"
    call GENQMD(neqnsz,IA,JA,IPErm,INVp,DEG,MARker,RCHset,NBRhd,QSIze,QLInk,NOFsub,ADJncy)
    do
      !   build up the parent vector parent vector will be saved in
      !   work2 for a while
      call GENPAQ(IA,JA,INVp,IPErm,MARker,NEQns,RCHset)
      !
      !   build up the binary tree
      !
      lbtree = lwk3 - 2*NEQns
      last = lbtree
      !rmem
      allocate (BTRee(2*(NEQns+1)),stat=IERror)
      if ( IERror/=0 ) stop "ALLOCATION ERROR, btree: SUB. matini"
      call GENBTQ(INVp,MARker,BTRee,ZPIv,izz,NEQns)
      !
      !   rotate the binary tree to avoid a zero pivot
      !
      if ( izz==0 ) then
        !
        !   post ordering
        !
        lpordr = last - neqns1
        last = lpordr
        !rmem
        allocate (PARent(NEQns),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, parent: SUB. matini.f"
        allocate (NCH(NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, nch: SUB. matini.f"
        allocate (PORdr(NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, pordr: SUB. matini.f"
        call POSORD(PARent,BTRee,INVp,IPErm,PORdr,NCH,NEQns,DEG,MARker,RCHset)
        !
        !   generate skelton graph
        !
        lleaf = last - NTTbr
        lxleaf = lleaf - neqns1
        ladp = lxleaf - neqns1
        last = ladp
        !rmem
        allocate (ADJncp(NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, adjncp: SUB. matini.f"
        allocate (XLEaf(NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, xleaf: SUB. matini.f"
        allocate (LEAf(NTTbr),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, leaf: SUB. matini.f"
        call GNLEAF(IA,JA,INVp,IPErm,NCH,ADJncp,XLEaf,LEAf,NEQns,lnleaf)
        call FORPAR(NEQns,PARent,NCH,NSTop)
        !*Deallocation of work arrays
        deallocate (IA)
        deallocate (JA)
        deallocate (DEG)
        deallocate (MARker)
        deallocate (RCHset)
        deallocate (NBRhd)
        deallocate (QSIze)
        deallocate (QLInk)
        deallocate (ADJncy)
        deallocate (ZPIv)
        !*Nullify pointers
        nullify (IA)
        nullify (JA)
        nullify (DEG)
        nullify (MARker)
        nullify (RCHset)
        nullify (NBRhd)
        nullify (QSIze)
        nullify (QLInk)
        nullify (ADJncy)
        nullify (ZPIv)
        !
        !   build up xlnzr,colno  (this is the symbolic fct.)
        !
        maxl = lxleaf - (left+neqns1)
        allocate (XLNzr(NEQns+1),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, xlnzr: SUB. matini.f"
        call PRE_GNCLNO(PARent,XLEaf,LEAf,XLNzr,NEQns,NSTop,lncol,ir1)
        allocate (COLno(lncol),stat=IERror)
        if ( IERror/=0 ) stop "ALLOCATION ERROR, colno: SUB. matini.f"
        call GNCLNO(PARent,XLEaf,LEAf,XLNzr,COLno,NEQns,NSTop,lncol,ir1)
        !*Deallocate work arrays
        deallocate (PORdr)
        deallocate (ADJncp)
        deallocate (XLEaf)
        deallocate (LEAf)
        deallocate (BTRee)
        !*Nullify pointers
        nullify (PORdr)
        nullify (ADJncp)
        nullify (XLEaf)
        nullify (LEAf)
        nullify (BTRee)
        !rmem
        left = (left+neqns1) + lncol
        !rmiv
        LEN_dsln = (NEQns-NSTop+1)*(NEQns-NSTop)/2

        !Scalar assignments
        LEN_colno = lncol
        !
        !   area for REAL(kind=8) values
        !
        if ( mod(left,2)==0 ) left = left + 1

        !rmiv
        TOTal = left
        !rmiv
        STAge = 10
        IALoc = 5*NEQns + lncol + 1
        exit
      else
        if ( izz0==0 ) izz0 = izz
        if ( izz0/=izz ) then
          call BRINGU(ZPIv,IPErm,INVp,MARker,izz,NEQns,IRR)
        else
          lwk4 = last - neqns1
          lwk5 = lwk4 - neqns1
          last = lwk5
          call ROTATE(IA,JA,INVp,IPErm,MARker,BTRee,izz,NEQns,NBRhd,QSIze,IRR)
        endif
      endif
    enddo
  end subroutine MATINI

  !======================================================================!
  !> @brief NUFORM
  !======================================================================!
  subroutine NUFORM(hecMESH,hecMAT,Ir)
    use HECMW_UTIL
    implicit none
    !------
    type (HECMWST_LOCAL_MESH), intent(in)::hecMESH
    type (HECMWST_MATRIX), intent(in)::hecMAT
    integer(kind=kint), intent(out):: Ir
    !------
    integer:: i
    integer:: idbg
    integer:: ierr
    integer:: j
    integer:: k
    integer:: kk
    integer:: ntotal
    integer(kind=kint):: numnp
    integer(kind=kint):: ndof
    integer(kind=kint):: ndof2

    idbg = 0
    numnp = hecMAT%NP
    ndof = hecMESH%N_DOF
    ntotal = numnp*ndof

    !*NUFACT variables
    NEQns = numnp
    NDEg = ndof
    NTTbr = hecMAT%NP + hecMAT%NPL
    !+hecMAT%NPU if unsymmetric
    ISYm = 0

    !*Allocations
    allocate (val(NDEg*NDEg),stat=ierr)
    if ( ierr/=0 ) stop "Allocation error:val"
    write (6,*) "nuform:stage = ", STAge
    kk = 0
    ndof2 = ndof*ndof

    do j = 1, numnp
      !*Diagonal
      kk = kk + 1
      call VLCPY(val,hecMAT%D(ndof2*(j-1)+1:ndof2*j),ndof)
      call STAIJ1(0,j,j,val,Ir)

      do i = 1, ndof
        if ( val((i-1)*ndof+i)<=0 ) write (idbg,*) 'j,j,val:', j, i, val((i-1)*ndof+i)
      enddo

      !*Lower
      do k = hecMAT%INDEXL(j-1) + 1, hecMAT%INDEXL(j)
        i = hecMAT%ITEML(k)
        kk = kk + 1
        call VLCPY(val,hecMAT%AL(ndof2*(k-1)+1:ndof2*k),ndof)
        call STAIJ1(0,j,i,val,Ir)
      enddo
    enddo

    deallocate (val)
  end subroutine NUFORM

  !======================================================================!
  !> @brief NUFCT0 performs Cholesky factorization
  !          if(iv(22).eq.0)    normal type
  !          if(iv(22).gt.0)    code generation type
  !======================================================================!
  subroutine NUFCT0(Ir)
    use HECMW_UTIL
    implicit none
    !------
    integer, intent(out):: Ir
    !------
    integer:: ndeg2
    integer:: ndegl

    if ( STAge/=20 ) then
      print *, '*********Setting Stage 40!*********'
      Ir = 40
      return
    else
      Ir = 0
    endif

    allocate (TEMp(NDEg*NDEg*NEQns),stat=IRR)
    if ( IRR/=0 ) then
      write (*,*) '##Error : Not enough memory'
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
      !stop
    endif
    allocate (INDx(NEQns),stat=IRR)
    if ( IRR/=0 ) then
      write (*,*) '##Error : Not enough memory'
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
      !stop
    endif
    !rmiv
    ndegl = NDEg*(NDEg+1)
    ndegl = ndegl/2
    ndeg2 = NDEg*NDEg
    !rmiv
    select case (NDEg)
      case (1)
        call NUFCT(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
      case (2)
        call NUFCT2(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
      case (3)
        call NUFCT3(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
      case (6)
        call NUFCT6(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
      case default
        call NUFCTX(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,NDEg,ndegl,Ir)
    endselect

    STAge = 30
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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Ir
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Zln(*)
    real(kind=8), intent(out):: Diag(*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Dsln(*)
    !------
    integer:: ic
    integer:: l
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5
    real(kind=8):: tt

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
      call sum(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
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
    call SUM3(Neqns-Nstop+1,Dsln,Diag(Nstop),Indx,Temp)
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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Ir
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Zln(4,*)
    real(kind=8), intent(out):: Diag(3,*)
    real(kind=8), intent(out):: Temp(4,*)
    real(kind=8), intent(out):: Dsln(4,*)
    !------
    integer:: ic
    integer:: l
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5
    real(kind=8):: tt

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Ir
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Zln(9,*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Dsln(9,*)
    !------
    integer:: ic
    integer:: l
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5
    real(kind=8):: tt

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Zln(36,*)
    real(kind=8), intent(out):: Diag(21,*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Dsln(36,*)
    !------
    integer:: ic
    integer:: Ir
    integer:: l
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5
    real(kind=8):: tt

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
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Zln(Ndeg*Ndeg,*)
    real(kind=8), intent(out):: Diag(Ndegl,*)
    real(kind=8), intent(out):: Temp(Ndeg*Ndeg,*)
    real(kind=8), intent(out):: Dsln(Ndeg*Ndeg,*)
    !------
    integer:: ic
    integer:: Ir
    integer:: l
    real(kind=8):: t1
    real(kind=8):: t2
    real(kind=8):: t3
    real(kind=8):: t4
    real(kind=8):: t5
    real(kind=8):: tt
    real(kind=8):: zz(100)
    real(kind=8):: t(100)

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
  !> @brief NUSOL0 performs forward elimination and backward substitution
  !     (i/o)
  !           r_h_s    on entry     right hand side vector
  !                    on exit      solution vector
  !           iv       communication array
  !======================================================================!
  subroutine NUSOL0(R_h_s,Ir)
    use HECMW_UTIL
    implicit none
    !------
    integer, intent(out):: Ir
    real(kind=8), intent(out):: R_h_s(*)
    !------
    integer:: lwk
    integer:: ndegl
    real(kind=8), pointer :: wk(:)

    if ( STAge/=30 .and. STAge/=40 ) then
      Ir = 50
      return
    else
      Ir = 0
    endif
    lwk = TOTal

    allocate (wk(NDEg*NEQns),stat=IERror)
    if ( IERror/=0 ) then
      write (*,*) "##Error: not enough memory"
      call HECMW_ABORT(HECMW_COMM_GET_COMM())
    endif
    !rmiv
    ndegl = NDEg*(NDEg+1)
    ndegl = ndegl/2

    select case( NDEg )
      case (1)
        call NUSOL1(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
      case (2)
        call NUSOL2(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
      case (3)
        call NUSOL3(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
      case (6)
        call NUSOLX(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop,NDEg,ndegl)
      case default
        call NUSOLX(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop,NDEg,ndegl)
    endselect

    STAge = 40
    deallocate (wk)
  end subroutine NUSOL0

  !======================================================================!
  !> @brief NUSOL1 performs forward elimination and backward substitution
  !======================================================================!
  subroutine NUSOL1(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Iperm(*)
    real(kind=8), intent(in):: Zln(*)
    real(kind=8), intent(in):: Diag(*)
    real(kind=8), intent(in):: Dsln(*)
    real(kind=8), intent(out):: B(*)
    real(kind=8), intent(out):: Wk(*)
    !------
    integer:: i
    integer:: j
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
        Wk(i) = Wk(i) - DDOT(Wk(Nstop),Dsln(joc),i-Nstop)
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
    ! permutaion
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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Iperm(*)
    real(kind=8), intent(in):: Zln(4,*)
    real(kind=8), intent(in):: Diag(3,*)
    real(kind=8), intent(in):: Dsln(4,*)
    real(kind=8), intent(out):: B(2,*)
    real(kind=8), intent(out):: Wk(2,*)
    !------
    integer:: i
    integer:: j
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
    ! permutaion
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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Iperm(*)
    real(kind=8), intent(in):: Zln(9,*)
    real(kind=8), intent(in):: Diag(6,*)
    real(kind=8), intent(in):: Dsln(9,*)
    real(kind=8), intent(out):: B(3,*)
    real(kind=8), intent(out):: Wk(3,*)
    !------
    integer:: i
    integer:: j
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
    ! permutaion
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
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Iperm(*)
    real(kind=8), intent(in):: Zln(Ndeg,Ndeg,*)
    real(kind=8), intent(in):: Diag(Ndegl,*)
    real(kind=8), intent(out):: B(Ndeg,*)
    real(kind=8), intent(out):: Wk(Ndeg,*)
    real(kind=8), intent(out):: Dsln(Ndeg,Ndeg,*)
    !------
    integer:: i
    integer:: j
    integer:: joc
    integer:: joc1
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    integer:: loc1
    integer:: locd
    integer:: m
    integer:: n

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
    ! permutaion
    do l = 1, Ndeg
      B(l,Iperm(1:Neqns)) = Wk(l,1:Neqns)
    enddo
  end subroutine NUSOLX

  !======================================================================!
  !> @brief ZPIVOT
  !======================================================================!
  subroutine ZPIVOT(Neqns,Neqnsz,Nttbr,Jcol,Irow,Zpiv,Ir)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Nttbr
    integer, intent(in):: Jcol(*)
    integer, intent(in):: Irow(*)
    integer, intent(out):: Ir
    integer, intent(out):: Neqnsz
    integer, intent(out):: Zpiv(*)
    !------
    integer:: i
    integer:: j
    integer:: l

    Ir = 0
    Zpiv(1:Neqns) = 1

    do l = 1, Nttbr
      i = Irow(l)
      j = Jcol(l)
      if ( i<=0 .or. j<=0 ) then
        Ir = -1
        goto 100
      elseif ( i>Neqns .or. j>Neqns ) then
        Ir = 1
        goto 100
      endif
      if ( i==j ) Zpiv(i) = 0
    enddo
    do i = Neqns, 1, -1
      if ( Zpiv(i)==0 ) then
        Neqnsz = i
        exit
      endif
    enddo
    100 if ( IDBg/=0 ) write (6,"(20I3)") (Zpiv(i),i=1,Neqns)
  end subroutine ZPIVOT

  !======================================================================!
  !> @brief STSMAT
  !======================================================================!
  subroutine STSMAT(Neqns,Nttbr,Irow,Jcol,Jcpt,Jcolno)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Nttbr
    integer, intent(in):: Irow(*)
    integer, intent(in):: Jcol(*)
    integer, intent(out):: Jcpt(*)
    integer, intent(out):: Jcolno(*)
    !------
    integer:: i
    integer:: j
    integer:: joc
    integer:: k
    integer:: l
    integer:: locr

    Jcpt(1:2*Nttbr) = 0
    Jcolno(1:2*Nttbr) = 0
    do i = 1, Neqns
      Jcpt(i) = i + Neqns
      Jcolno(i+Neqns) = i
    enddo

    k = 2*Neqns

    do l = 1, Nttbr
      i = Irow(l)
      j = Jcol(l)
      if ( i/=j ) then
        joc = Jcpt(i)
        locr = i
        do while ( joc/=0 )
          if ( Jcolno(joc)==j ) goto 100
          if ( Jcolno(joc)>j ) then
            k = k + 1
            Jcpt(locr) = k
            Jcpt(k) = joc
            Jcolno(k) = j
            goto 20
          endif
          locr = joc
          joc = Jcpt(joc)
        enddo
        k = k + 1
        Jcpt(locr) = k
        Jcolno(k) = j

        20      joc = Jcpt(j)
        locr = j
        do while ( joc/=0 )
          if ( Jcolno(joc)==i ) goto 100
          if ( Jcolno(joc)>i ) then
            k = k + 1
            Jcpt(locr) = k
            Jcpt(k) = joc
            Jcolno(k) = i
            goto 100
          endif
          locr = joc
          joc = Jcpt(joc)
        enddo
        k = k + 1
        Jcpt(locr) = k
        Jcolno(k) = i
      endif
      100 enddo
      if ( IDBg/=0 ) then
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
    integer, intent(in):: Neqns
    integer, intent(in):: Jcpt(*)
    integer, intent(in):: Jcolno(*)
    integer, intent(out):: Ia(*)
    integer, intent(out):: Ja(*)
    !------
    integer:: i
    integer:: ii
    integer:: joc
    integer:: k
    integer:: l

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
    if ( IDBg/=0 ) then
      write (6,*) 'ia '
      write (6,"(10I7)") (Ia(i),i=1,Neqns)
      write (6,*) 'ja '
      write (6,"(10I7)") (Ja(i),i=1,Ia(Neqns+1))
    endif
  end subroutine STIAJA

  !======================================================================!
  !> @brief GENQMD
  !======================================================================!
  subroutine GENQMD(Neqns,Xadj,Adj0,Perm,Invp,Deg,Marker,Rchset,Nbrhd,Qsize,Qlink,Nofsub,Adjncy)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Adj0(*)
    integer, intent(in):: Xadj(*)
    integer, intent(out):: Nofsub
    integer, intent(out):: Adjncy(*)
    integer, intent(out):: Perm(*)
    integer, intent(out):: Invp(*)
    integer, intent(out):: Deg(*)
    integer, intent(out):: Marker(*)
    integer, intent(out):: Rchset(*)
    integer, intent(out):: Nbrhd(*)
    integer, intent(out):: Qsize(*)
    integer, intent(out):: Qlink(*)
    !------
    integer:: inode
    integer:: ip
    integer:: irch
    integer:: j
    integer:: mindeg
    integer:: ndeg
    integer:: nhdsze
    integer:: node
    integer:: np
    integer:: num
    integer:: nump1
    integer:: nxnode
    integer:: rchsze
    integer:: search
    integer:: thresh

    mindeg = Neqns
    Nofsub = 0
    Adjncy(1:Xadj(Neqns+1) - 1) = Adj0(1:Xadj(Neqns+1) - 1)
    do node = 1, Neqns
      Perm(node) = node
      Invp(node) = node
      Marker(node) = 0
      Qsize(node) = 1
      Qlink(node) = 0
      ndeg = Xadj(node+1) - Xadj(node)
      Deg(node) = ndeg
      if ( ndeg<mindeg ) mindeg = ndeg
    enddo

    num = 0
    100 search = 1
    thresh = mindeg
    mindeg = Neqns
    200 nump1 = num + 1
    if ( nump1>search ) search = nump1
    do j = search, Neqns
      node = Perm(j)
      if ( Marker(node)>=0 ) then
        ndeg = Deg(node)
        if ( ndeg<=thresh ) goto 300
        if ( ndeg<mindeg ) mindeg = ndeg
      endif
    enddo
    goto 100

    300 search = j
    Nofsub = Nofsub + Deg(node)
    Marker(node) = 1
    call QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
    nxnode = node
    do
      num = num + 1
      np = Invp(nxnode)
      ip = Perm(num)
      Perm(np) = ip
      Invp(ip) = np
      Perm(num) = nxnode
      Invp(nxnode) = num
      Deg(nxnode) = -1
      nxnode = Qlink(nxnode)
      if ( nxnode<=0 ) then
        if ( rchsze>0 ) then
          call QMDUPD(Xadj,Adjncy,rchsze,Rchset,Deg,Qsize,Qlink,Marker,Rchset(rchsze+1),Nbrhd(nhdsze+1))
          Marker(node) = 0
          do irch = 1, rchsze
            inode = Rchset(irch)
            if ( Marker(inode)>=0 ) then
              Marker(inode) = 0
              ndeg = Deg(inode)
              if ( ndeg<mindeg ) mindeg = ndeg
              if ( ndeg<=thresh ) then
                mindeg = thresh
                thresh = ndeg
                search = Invp(inode)
              endif
            endif
          enddo
          if ( nhdsze>0 ) call QMDOT(node,Xadj,Adjncy,Marker,rchsze,Rchset,Nbrhd)
        endif
        if ( num>=Neqns ) exit
        goto 200
      endif
    enddo
  end subroutine GENQMD

  !======================================================================!
  !> @brief GENPAQ
  !======================================================================!
  subroutine GENPAQ(Xadj,Adjncy,Invp,Iperm,Parent,Neqns,Ancstr)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Xadj(*)
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: Invp(*)
    integer, intent(in):: Iperm(*)
    integer, intent(out):: Parent(*)
    integer, intent(out):: Ancstr(*)
    !------
    integer:: i
    integer:: ip
    integer:: it
    integer:: k
    integer:: l

    do i = 1, Neqns
      Parent(i) = 0
      Ancstr(i) = 0
      ip = Iperm(i)
      do k = Xadj(ip), Xadj(ip+1) - 1
        l = Invp(Adjncy(k))
        if ( l<i ) then
          do while ( Ancstr(l)/=0 )
            if ( Ancstr(l)==i ) goto 50
            it = Ancstr(l)
            Ancstr(l) = i
            l = it
          enddo
          Ancstr(l) = i
          Parent(l) = i
        endif
        50    enddo
      enddo
      do i = 1, Neqns
        if ( Parent(i)==0 ) Parent(i) = Neqns + 1
      enddo
      Parent(Neqns+1) = 0

      if ( IDBg1/=0 ) then
        write (6,"(' parent')")
        write (6,"(2I6)") (i,Parent(i),i=1,Neqns)
      endif
  end subroutine GENPAQ

  !======================================================================!
  !> @brief GENBTQ
  !======================================================================!
  subroutine GENBTQ(Invp,Parent,Btree,Zpiv,Izz,Neqns)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Parent(*)
    integer, intent(in):: Invp(*)
    integer, intent(in):: Zpiv(*)
    integer, intent(out):: Btree(2,*)
    !------
    integer:: i
    integer:: ib
    integer:: inext
    integer:: ip
    integer:: Izz

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
    do i = 1, Neqns
      if ( Zpiv(i)/=0 ) then
        if ( Btree(1,Invp(i))==0 ) then
          Izz = i
          goto 100
        endif
      endif
    enddo
    Izz = 0

    100 continue
    if ( IDBg1/=0 ) then
      write (6,"(' binary tree')")
      write (6,"(i6,'(',2I6,')')") (i,Btree(1,i),Btree(2,i),i=1,Neqns)
      write (6,"(' the first zero pivot is ',i4)") Izz
    endif
  end subroutine GENBTQ

  !======================================================================!
  !> @brief POSORD
  !======================================================================!
  subroutine POSORD(Parent,Btree,Invp,Iperm,Pordr,Nch,Neqns,Iw,Qarent,Mch)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Btree(2,*)
    integer, intent(in):: Qarent(*)
    integer, intent(out):: Parent(*)
    integer, intent(out):: Pordr(*)
    integer, intent(out):: Nch(*)
    integer, intent(out):: Invp(*)
    integer, intent(out):: Iperm(*)
    integer, intent(out):: Iw(*)
    integer, intent(out):: Mch(0:Neqns+1)
    !------
    integer:: i
    integer:: ii
    integer:: invpos
    integer:: ipinv
    integer:: joc
    integer:: l
    integer:: locc
    integer:: locp

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
            if ( IDBg1/=0 ) then
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
  end subroutine POSORD

  !======================================================================!
  !> @brief GNLEAF
  !======================================================================!
  subroutine GNLEAF(Xadj,Adjncy,Invp,Iperm,Nch,Adjncp,Xleaf,Leaf,Neqns,Lnleaf)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Xadj(*)
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: Nch(*)
    integer, intent(in):: Invp(*)
    integer, intent(in):: Iperm(*)
    integer, intent(out):: Lnleaf
    integer, intent(out):: Adjncp(*)
    integer, intent(out):: Xleaf(*)
    integer, intent(out):: Leaf(*)
    !------
    integer:: i
    integer:: ik
    integer:: ip
    integer:: iq
    integer:: istart
    integer:: k
    integer:: l
    integer:: lc
    integer:: lc1
    integer:: m

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
        call QQSORT(Adjncp(istart+1),m)
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

    if ( IDBg1/=0 ) then
      write (6,"(' xleaf')")
      write (6,"(10I6)") (Xleaf(i),i=1,Neqns+1)
      write (6,"(' leaf (len = ',i6,')')") Lnleaf
      write (6,"(10I6)") (Leaf(i),i=1,Lnleaf)
    endif

    return
  end subroutine GNLEAF

  !======================================================================!
  !> @brief FORPAR
  !======================================================================!
  subroutine FORPAR(Neqns,Parent,Nch,Nstop)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Parent(*)
    integer, intent(out):: Nstop
    integer, intent(out):: Nch(*)
    !------
    integer:: i
    integer:: idens
    integer:: ii

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Parent(*)
    integer, intent(in):: Xleaf(*)
    integer, intent(in):: Leaf(*)
    integer, intent(out):: Lncol
    integer, intent(out):: Xlnzr(*)
    !------
    integer:: i
    integer:: Ir
    integer:: j
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    integer:: nc
    integer:: nxleaf

    nc = 0
    Ir = 0
    l = 1
    do i = 1, Neqns
      Xlnzr(i) = l
      ks = Xleaf(i)
      ke = Xleaf(i+1) - 1
      if ( ke>=ks ) then
        nxleaf = Leaf(ks)
        do k = ks, ke - 1
          j = nxleaf
          nxleaf = Leaf(k+1)
          do while ( j<nxleaf )
            if ( j>=Nstop ) goto 100
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
      100 enddo
      Xlnzr(Neqns+1) = l
      Lncol = l - 1
  end subroutine PRE_GNCLNO

  !======================================================================!
  !> @brief GNCLNO
  !======================================================================!
  subroutine GNCLNO(Parent,Xleaf,Leaf,Xlnzr,Colno,Neqns,Nstop,Lncol,Ir)
    implicit none
    !------
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Parent(*)
    integer, intent(in):: Xleaf(*)
    integer, intent(in):: Leaf(*)
    integer, intent(out):: Ir
    integer, intent(out):: Lncol
    integer, intent(out):: Xlnzr(*)
    integer, intent(out):: Colno(*)
    !------
    integer:: i
    integer:: j
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    integer:: nc
    integer:: nxleaf

    nc = 0
    Ir = 0
    l = 1
    do i = 1, Neqns
      Xlnzr(i) = l
      ks = Xleaf(i)
      ke = Xleaf(i+1) - 1
      if ( ke>=ks ) then
        nxleaf = Leaf(ks)
        do k = ks, ke - 1
          j = nxleaf
          nxleaf = Leaf(k+1)
          do while ( j<nxleaf )
            if ( j>=Nstop ) goto 100
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
      100 enddo
      Xlnzr(Neqns+1) = l
      Lncol = l - 1

      if ( IDBg1/=0 ) then
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
    integer, intent(in):: Izz
    integer, intent(in):: Neqns
    integer, intent(in):: Zpiv(*)
    integer, intent(in):: Parent(*)
    integer, intent(out):: Irr
    integer, intent(out):: Iperm(*)
    integer, intent(out):: Invp(*)
    !------
    integer:: i
    integer:: ib
    integer:: ib0
    integer:: ibp
    integer:: idbg
    integer:: izzp

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
  subroutine ROTATE(Xadj,Adjncy,Invp,Iperm,Parent,Btree,Izz,Neqns,Anc,Adjt,Irr)
    implicit none
    !------
    integer, intent(in):: Izz
    integer, intent(in):: Neqns
    integer, intent(in):: Xadj(*)
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: Parent(*)
    integer, intent(in):: Btree(2,*)
    integer, intent(out):: Irr
    integer, intent(out):: Invp(*)
    integer, intent(out):: Iperm(*)
    integer, intent(out):: Anc(*)
    integer, intent(out):: Adjt(*)
    !------
    integer:: i
    integer:: iy
    integer:: izzz
    integer:: joc
    integer:: k
    integer:: l
    integer:: ll
    integer:: locc
    integer:: nanc

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

    100 continue
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
                goto 100
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
              do
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
                      goto 105
                    else
                      locc = Btree(2,joc)
                      if ( locc/=0 ) exit
                      joc = Parent(joc)
                    endif
                  enddo
                endif
              enddo
            endif
            !
            ! set iperm
            !
            105         do i = 1, Neqns
            Iperm(Invp(i)) = i
        enddo

        if ( IDBg1/=0 ) write (6,"(10I6)") (Invp(i),i=1,Neqns)
        return
      else
        locc = Btree(2,joc)
        if ( locc/=0 ) exit
        joc = Parent(joc)
      endif
    enddo
  endif
enddo
  end subroutine ROTATE

  !======================================================================!
  !> @brief ADDR0
  !======================================================================!
  subroutine ADDR0(Isw,I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ndeg2,Ndeg2l,Ir)
    implicit none
    !------
    integer, intent(in):: I
    integer, intent(in):: Isw
    integer, intent(in):: J
    integer, intent(in):: Ndeg2
    integer, intent(in):: Ndeg2l
    integer, intent(in):: Nstop
    integer, intent(in):: Invp(*)
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Aij(Ndeg2)
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Zln(Ndeg2,*)
    real(kind=8), intent(out):: Diag(Ndeg2l,*)
    real(kind=8), intent(out):: Dsln(Ndeg2,*)
    !------
    integer:: i0
    integer:: ii
    integer:: itrans
    integer:: j0
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer, parameter:: idbg = 0

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
    integer, intent(in):: I
    integer, intent(in):: J
    integer, intent(in):: Nstop
    integer, intent(in):: Invp(*)
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Aij(9)
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Zln(9,*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Dsln(9,*)
    !------
    integer:: i0
    integer:: ii
    integer:: itrans
    integer:: j0
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    integer, parameter:: idbg = 0
    integer, parameter:: ndeg2 = 9
    integer, parameter:: ndeg2l = 6

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
    integer, intent(in):: I
    integer, intent(in):: J
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndeg2l
    integer, intent(in):: Nstop
    integer, intent(in):: Invp(*)
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Aij(Ndeg,Ndeg)
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Zln(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Diag(Ndeg2l,*)
    real(kind=8), intent(out):: Dsln(Ndeg,Ndeg,*)
    !------
    integer:: i0
    integer:: ii
    integer:: itrans
    integer:: j0
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    integer:: m
    integer:: n
    integer, parameter:: idbg = 0

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
  !> @brief D2DOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D2DOT(T,A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(4,*)
    real(kind=8), intent(in):: B(4,*)
    real(kind=8), intent(out):: T(4)
    !------
    integer:: jj

    T(1:4) = 0.0D0

    do jj = 1, N
      T(1) = T(1) + A(1,jj)*B(1,jj) + A(3,jj)*B(3,jj)
      T(2) = T(2) + A(2,jj)*B(1,jj) + A(4,jj)*B(3,jj)
      T(3) = T(3) + A(1,jj)*B(2,jj) + A(3,jj)*B(4,jj)
      T(4) = T(4) + A(2,jj)*B(2,jj) + A(4,jj)*B(4,jj)
    enddo
  end subroutine D2DOT

  !======================================================================!
  !> @brief D2SDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D2SDOT(Wi,A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(2,*)
    real(kind=8), intent(in):: B(4,*)
    real(kind=8), intent(out):: Wi(2)
    !------
    integer:: jj

    do jj = 1, N
      Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(3,jj)
      Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(4,jj)
    enddo
  end subroutine D2SDOT

  !======================================================================!
  !> @brief D3DOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D3DOT(T,A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(9,*)
    real(kind=8), intent(in):: B(9,*)
    real(kind=8), intent(out):: T(9)
    !------
    integer:: jj

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
    integer, intent(in):: N
    real(kind=8), intent(in):: A(9,*)
    real(kind=8), intent(in):: B(9,*)
    real(kind=8), intent(out):: T(6)
    !------
    integer:: jj

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
  !> @brief D3SDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine D3SDOT(Wi,A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(3,*)
    real(kind=8), intent(in):: B(9,*)
    real(kind=8), intent(out):: Wi(3)
    !------
    integer:: jj

    do jj = 1, N
      Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(4,jj) - A(3,jj)*B(7,jj)
      Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(5,jj) - A(3,jj)*B(8,jj)
      Wi(3) = Wi(3) - A(1,jj)*B(3,jj) - A(2,jj)*B(6,jj) - A(3,jj)*B(9,jj)
    enddo
  end subroutine D3SDOT

  !======================================================================!
  !> @brief DDOT performs inner product of sparse vectors
  !======================================================================!
  real(kind=8) function DDOT(A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(N)
    real(kind=8), intent(in):: B(N)
    !------
    integer:: i
    real(kind=8):: s

    s = 0.0D0
    do i = 1, N
      s = s + A(i)*B(i)
    enddo
    DDOT = s
  end function DDOT

  !======================================================================!
  !> @brief DXDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine DXDOT(Ndeg,T,A,B,L)
    implicit none
    !------
    integer, intent(in):: L
    integer, intent(in):: Ndeg
    real(kind=8), intent(in):: A(Ndeg,Ndeg,*)
    real(kind=8), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: T(Ndeg,Ndeg)
    !------
    integer:: jj
    integer:: k
    integer:: m
    integer:: n

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
    integer, intent(in):: L
    integer, intent(in):: Ndeg
    real(kind=8), intent(in):: A(Ndeg,Ndeg,*)
    real(kind=8), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: T(Ndeg,Ndeg)
    !------
    integer:: jj
    integer:: k
    integer:: m
    integer:: n

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
  !> @brief DXSDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine DXSDOT(Ndeg,Wi,A,B,N)
    implicit none
    !------
    integer, intent(in):: Ndeg
    real(kind=8), intent(in):: A(Ndeg,*)
    real(kind=8), intent(in):: B(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Wi(Ndeg)
    integer, intent(inout):: N
    !------
    integer:: jj
    integer:: m

    do jj = 1, N
      do m = 1, Ndeg
        do N = 1, Ndeg
          Wi(N) = Wi(N) - B(N,m,jj)*A(m,jj)
        enddo
      enddo
    enddo
  end subroutine DXSDOT

  !======================================================================!
  !> @brief S3PDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine S3PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer, intent(in):: Ke
    integer, intent(in):: Ks
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Zln(9,*)
    real(kind=8), intent(in):: B(3,*)
    real(kind=8), intent(out):: Bi(3)
    !------
    integer:: j
    integer:: jj

    do jj = Ks, Ke
      j = Colno(jj)
      Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(4,jj)*B(2,j) - Zln(7,jj)*B(3,j)
      Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(5,jj)*B(2,j) - Zln(8,jj)*B(3,j)
      Bi(3) = Bi(3) - Zln(3,jj)*B(1,j) - Zln(6,jj)*B(2,j) - Zln(9,jj)*B(3,j)
    enddo
  end subroutine S3PDOT

  !======================================================================!
  !> @brief S2PDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine S2PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer, intent(in):: Ke
    integer, intent(in):: Ks
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Zln(4,*)
    real(kind=8), intent(in):: B(2,*)
    real(kind=8), intent(out):: Bi(2)
    !------
    integer:: j
    integer:: jj

    do jj = Ks, Ke
      j = Colno(jj)
      Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(3,jj)*B(2,j)
      Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(4,jj)*B(2,j)
    enddo
  end subroutine S2PDOT

  !======================================================================!
  !> @brief S6PDOT  performs inner product of sparse vectors
  !======================================================================!
  subroutine S6PDOT(Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer, intent(in):: Ke
    integer, intent(in):: Ks
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Zln(36,*)
    real(kind=8), intent(in):: B(6,*)
    real(kind=8), intent(out):: Bi(6)
    !------
    integer:: j
    integer:: jj

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
  !> @brief SPDOT2 performs inner product of sparse vectors
  !======================================================================!
  real(kind=8) function SPDOT2(B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer, intent(in):: Ke
    integer, intent(in):: Ks
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Zln(*)
    real(kind=8), intent(in):: B(*)
    !------
    integer:: j
    integer:: jj
    real(kind=8):: s

    s = 0.0D0
    do jj = Ks, Ke
      j = Colno(jj)
      s = s + Zln(jj)*B(j)
    enddo
    SPDOT2 = s
  end function SPDOT2

  !======================================================================!
  !> @brief SXPDOT performs inner product of sparse vectors
  !======================================================================!
  subroutine SXPDOT(Ndeg,Bi,B,Zln,Colno,Ks,Ke)
    implicit none
    !------
    integer, intent(in):: Ke
    integer, intent(in):: Ks
    integer, intent(in):: Ndeg
    integer, intent(in):: Colno(*)
    real(kind=8), intent(in):: Zln(Ndeg,Ndeg,*)
    real(kind=8), intent(in):: B(Ndeg,*)
    real(kind=8), intent(out):: Bi(Ndeg)
    !------
    integer:: j
    integer:: jj
    integer:: m
    integer:: n

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
  !> @brief INV2
  !======================================================================!
  subroutine INV2(Dsln,Ir)
    implicit none
    !------
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Dsln(3)
    !------
    real(kind=8):: t

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
    real(kind=8), intent(in):: Diag(3)
    real(kind=8), intent(in):: Zz(4)
    real(kind=8), intent(out):: Zln(4)
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
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Dsln(6)
    !------
    real(kind=8):: t(2)

    Ir = 0
    if ( dabs(Dsln(1))<RMIn ) goto 100
    Dsln(1) = 1.0D0/Dsln(1)
    t(1) = Dsln(2)*Dsln(1)
    Dsln(3) = Dsln(3) - t(1)*Dsln(2)
    Dsln(2) = t(1)
    if ( dabs(Dsln(3))<RMIn ) goto 100
    Dsln(3) = 1.0D0/Dsln(3)
    t(1) = Dsln(4)*Dsln(1)
    Dsln(5) = Dsln(5) - Dsln(2)*Dsln(4)
    t(2) = Dsln(5)*Dsln(3)
    Dsln(6) = Dsln(6) - t(1)*Dsln(4) - t(2)*Dsln(5)
    Dsln(4) = t(1)
    Dsln(5) = t(2)
    if ( dabs(Dsln(6))<RMIn ) goto 100
    Dsln(6) = 1.0D0/Dsln(6)
    return

    100 Dsln(1) = 1.0D0
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
    real(kind=8), intent(in):: Diag(6)
    real(kind=8), intent(in):: Zz(9)
    real(kind=8), intent(out):: Zln(9)
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
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Dsln(21)
    !------
    real(kind=8):: t(5)

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
    real(kind=8), intent(in):: Diag(21)
    real(kind=8), intent(in):: Zz(36)
    real(kind=8), intent(out):: Zln(36)
    !------
    integer:: i

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
    integer, intent(in):: Ndeg
    integer, intent(out):: Ir
    real(kind=8), intent(out):: Dsln(*)
    !------
    integer:: i
    integer:: j
    integer:: k
    integer:: k0
    integer:: l
    integer:: l0
    integer:: ld
    integer:: ll
    real(kind=8):: t
    real(kind=8):: tem

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
    integer, intent(in):: Ndeg
    real(kind=8), intent(in):: Diag(*)
    real(kind=8), intent(in):: Zz(Ndeg,Ndeg)
    real(kind=8), intent(out):: Zln(Ndeg,Ndeg)
    !------
    integer:: joc
    integer:: l
    integer:: loc1
    integer:: m
    integer:: n

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
  !> @brief QMDMRG
  !======================================================================!
  subroutine QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,Deg0,Nhdsze,Nbrhd,Rchset,Ovrlp)
    implicit none
    !------
    integer, intent(in):: Deg0
    integer, intent(in):: Nhdsze
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: Nbrhd(*)
    integer, intent(in):: Xadj(*)
    integer, intent(out):: Deg(*)
    integer, intent(out):: Qsize(*)
    integer, intent(out):: Qlink(*)
    integer, intent(out):: Marker(*)
    integer, intent(out):: Rchset(*)
    integer, intent(out):: Ovrlp(*)
    !------
    integer:: deg1
    integer:: head
    integer:: inhd
    integer:: iov
    integer:: irch
    integer:: j
    integer:: jstrt
    integer:: jstop
    integer:: link
    integer:: lnode
    integer:: mark
    integer:: mrgsze
    integer:: nabor
    integer:: node
    integer:: novrlp
    integer:: rchsze
    integer:: root

    if ( Nhdsze<=0 ) return
    do inhd = 1, Nhdsze
      root = Nbrhd(inhd)
      Marker(root) = 0
    enddo
    do inhd = 1, Nhdsze
      root = Nbrhd(inhd)
      Marker(root) = -1
      rchsze = 0
      novrlp = 0
      deg1 = 0
      do
        jstrt = Xadj(root)
        jstop = Xadj(root+1) - 1
        do j = jstrt, jstop
          nabor = Adjncy(j)
          root = -nabor
          if ( nabor<0 ) goto 50
          if ( nabor==0 ) exit
          mark = Marker(nabor)

          if ( mark>=0 ) then
            if ( mark<=0 ) then
              rchsze = rchsze + 1
              Rchset(rchsze) = nabor
              deg1 = deg1 + Qsize(nabor)
              Marker(nabor) = 1
            elseif ( mark<=1 ) then
              novrlp = novrlp + 1
              Ovrlp(novrlp) = nabor
              Marker(nabor) = 2
            endif
          endif
        enddo
        exit
        50    enddo
        head = 0
        mrgsze = 0
        do iov = 1, novrlp
          node = Ovrlp(iov)
          jstrt = Xadj(node)
          jstop = Xadj(node+1) - 1
          do j = jstrt, jstop
            nabor = Adjncy(j)
            if ( Marker(nabor)==0 ) then
              Marker(node) = 1
              goto 100
            endif
          enddo
          mrgsze = mrgsze + Qsize(node)
          Marker(node) = -1
          lnode = node
          do
            link = Qlink(lnode)
            if ( link<=0 ) then
              Qlink(lnode) = head
              head = node
              exit
            else
              lnode = link
            endif
          enddo
          100   enddo
          if ( head>0 ) then
            Qsize(head) = mrgsze
            Deg(head) = Deg0 + deg1 - 1
            Marker(head) = 2
          endif
          root = Nbrhd(inhd)
          Marker(root) = 0
          if ( rchsze>0 ) then
            do irch = 1, rchsze
              node = Rchset(irch)
              Marker(node) = 0
            enddo
          endif
        enddo
  end subroutine QMDMRG

  !======================================================================!
  !> @brief QMDOT
  !======================================================================!
  subroutine QMDOT(Root,Xadj,Adjncy,Marker,Rchsze,Rchset,Nbrhd)
    implicit none
    !------
    integer, intent(in):: Rchsze
    integer, intent(in):: Root
    integer, intent(in):: Marker(*)
    integer, intent(in):: Rchset(*)
    integer, intent(in):: Nbrhd(*)
    integer, intent(in):: Xadj(*)
    integer, intent(out):: Adjncy(*)
    !------
    integer:: inhd
    integer:: irch
    integer:: j
    integer:: jstrt
    integer:: jstop
    integer:: link
    integer:: nabor
    integer:: node

    irch = 0
    inhd = 0
    node = Root
    100 jstrt = Xadj(node)
    jstop = Xadj(node+1) - 2
    if ( jstop>=jstrt ) then
      do j = jstrt, jstop
        irch = irch + 1
        Adjncy(j) = Rchset(irch)
        if ( irch>=Rchsze ) goto 200
      enddo
    endif
    link = Adjncy(jstop+1)
    node = -link
    if ( link>=0 ) then
      inhd = inhd + 1
      node = Nbrhd(inhd)
      Adjncy(jstop+1) = -node
    endif
    goto 100

    200 Adjncy(j+1) = 0
    do irch = 1, Rchsze
      node = Rchset(irch)
      if ( Marker(node)>=0 ) then
        jstrt = Xadj(node)
        jstop = Xadj(node+1) - 1
        do j = jstrt, jstop
          nabor = Adjncy(j)
          if ( Marker(nabor)<0 ) then
            Adjncy(j) = Root
            exit
          endif
        enddo
      endif
    enddo
  end subroutine QMDOT

  !======================================================================!
  !> @brief QMDRCH
  !======================================================================!
  subroutine QMDRCH(Root,Xadj,Adjncy,Deg,Marker,Rchsze,Rchset,Nhdsze,Nbrhd)
    implicit none
    !------
    integer, intent(in):: Root
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: Deg(*)
    integer, intent(in):: Xadj(*)
    integer, intent(out):: Nhdsze
    integer, intent(out):: Rchsze
    integer, intent(out):: Marker(*)
    integer, intent(out):: Rchset(*)
    integer, intent(out):: Nbrhd(*)
    !------
    integer:: i
    integer:: istrt
    integer:: istop
    integer:: j
    integer:: jstrt
    integer:: jstop
    integer:: nabor
    integer:: node

    Nhdsze = 0
    Rchsze = 0
    istrt = Xadj(Root)
    istop = Xadj(Root+1) - 1
    if ( istop<istrt ) return
    do i = istrt, istop
      nabor = Adjncy(i)
      if ( nabor==0 ) return
      if ( Marker(nabor)==0 ) then
        if ( Deg(nabor)<0 ) then
          Marker(nabor) = -1
          Nhdsze = Nhdsze + 1
          Nbrhd(Nhdsze) = nabor
          do
            jstrt = Xadj(nabor)
            jstop = Xadj(nabor+1) - 1
            do j = jstrt, jstop
              node = Adjncy(j)
              nabor = -node
              if ( node<0 ) goto 10
              if ( node==0 ) exit
              if ( Marker(node)==0 ) then
                Rchsze = Rchsze + 1
                Rchset(Rchsze) = node
                Marker(node) = 1
              endif
            enddo
            exit
            10        enddo
          else
            Rchsze = Rchsze + 1
            Rchset(Rchsze) = nabor
            Marker(nabor) = 1
        endif
      endif
    enddo
  end subroutine QMDRCH

  !======================================================================!
  !> @brief QMDUPD
  !======================================================================!
  subroutine QMDUPD(Xadj,Adjncy,Nlist,List,Deg,Qsize,Qlink,Marker,Rchset,Nbrhd)
    implicit none
    !------
    integer, intent(in):: Nlist
    integer, intent(in):: Adjncy(*)
    integer, intent(in):: List(*)
    integer, intent(in):: Xadj(*)
    integer, intent(out):: Deg(*)
    integer, intent(out):: Marker(*)
    integer, intent(out):: Rchset(*)
    integer, intent(out):: Nbrhd(*)
    integer, intent(out):: Qsize(*)
    integer, intent(out):: Qlink(*)
    !------
    integer:: deg0
    integer:: deg1
    integer:: il
    integer:: inhd
    integer:: inode
    integer:: irch
    integer:: j
    integer:: jstrt
    integer:: jstop
    integer:: mark
    integer:: nabor
    integer:: nhdsze
    integer:: node
    integer:: rchsze

    if ( Nlist<=0 ) return
    deg0 = 0
    nhdsze = 0
    do il = 1, Nlist
      node = List(il)
      deg0 = deg0 + Qsize(node)
      jstrt = Xadj(node)
      jstop = Xadj(node+1) - 1
      do j = jstrt, jstop
        nabor = Adjncy(j)
        if ( Marker(nabor)==0 .and. Deg(nabor)<0 ) then
          Marker(nabor) = -1
          nhdsze = nhdsze + 1
          Nbrhd(nhdsze) = nabor
        endif
      enddo
    enddo

    if ( nhdsze>0 ) call QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,deg0,nhdsze,Nbrhd,Rchset,Nbrhd(nhdsze+1))
    do il = 1, Nlist
      node = List(il)
      mark = Marker(node)
      if ( mark<=1 .and. mark>=0 ) then
        call QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
        deg1 = deg0
        if ( rchsze>0 ) then
          do irch = 1, rchsze
            inode = Rchset(irch)
            deg1 = deg1 + Qsize(inode)
            Marker(inode) = 0
          enddo
        endif
        Deg(node) = deg1 - 1
        if ( nhdsze>0 ) then
          do inhd = 1, nhdsze
            inode = Nbrhd(inhd)
            Marker(inode) = 0
          enddo
        endif
      endif
    enddo
  end subroutine QMDUPD

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
    integer, intent(in):: Ik
    integer, intent(out):: Iw(*)
    !------
    integer:: itemp
    integer:: l
    integer:: m

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
  !> @brief S2UM
  !======================================================================!
  subroutine S2UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    implicit none
    !------
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Par(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Diag(3,*)
    real(kind=8), intent(out):: Temp(4,*)
    real(kind=8), intent(out):: Zln(4,*)
    !------
    integer:: ir
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    real(kind=8):: s(4)
    real(kind=8):: t(3)
    real(kind=8):: zz(4)

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Temp(4,*)
    real(kind=8), intent(out):: Zln(4,*)
    !------
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    real(kind=8):: s(4)

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(3,*)
    real(kind=8), intent(out):: Dsln(4,*)
    real(kind=8), intent(out):: Temp(4,*)
    real(kind=8), intent(out):: Zln(4,*)
    !------
    integer:: ic
    integer:: j
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
    integer, intent(in):: N
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(3,*)
    real(kind=8), intent(out):: Dsln(4,*)
    real(kind=8), intent(out):: Temp(4,*)
    !------
    integer:: i
    integer:: ir
    integer:: j
    integer:: joc
    real(kind=8):: t(4)

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Par(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Temp(9,*)
    real(kind=8), intent(out):: Zln(9,*)
    !------
    integer:: ir
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    integer:: l
    real(kind=8):: t(6)
    real(kind=8):: zz(9)

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Temp(9,*)
    real(kind=8), intent(out):: Zln(9,*)
    !------
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    real(kind=8):: s(9)

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Dsln(9,*)
    real(kind=8), intent(out):: Temp(Neqns,9)
    real(kind=8), intent(out):: Zln(9,*)
    !------
    integer:: ic
    integer:: j
    integer:: j1
    integer:: j2
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
    integer, intent(in):: N
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Dsln(9,*)
    real(kind=8), intent(out):: Temp(9,*)
    !------
    integer:: i
    integer:: ir
    integer:: j
    integer:: joc
    real(kind=8):: t(9)

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Par(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Diag(21,*)
    real(kind=8), intent(out):: Temp(36,*)
    real(kind=8), intent(out):: Zln(36,*)
    !------
    integer:: ir
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    integer:: l
    real(kind=8):: t(21)
    real(kind=8):: zz(36)

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Temp(9,*)
    real(kind=8), intent(out):: Zln(9,*)
    !------
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    integer:: l
    real(kind=8):: s(9)

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(21,*)
    real(kind=8), intent(out):: Dsln(36,*)
    real(kind=8), intent(out):: Temp(36,Neqns)
    real(kind=8), intent(out):: Zln(36,*)
    !------
    integer:: ic
    integer:: j
    integer:: j1
    integer:: j2
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks

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
    integer, intent(in):: N
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(6,*)
    real(kind=8), intent(out):: Dsln(9,*)
    real(kind=8), intent(out):: Temp(9,*)
    !------
    integer:: i
    integer:: ir
    integer:: j
    integer:: joc
    integer:: l
    real(kind=8):: t(9)

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
  !======================================================================!
  subroutine STAIJ1(Isw,I,J,Aij,Ir)
    implicit none
    !------
    integer, intent(in):: I
    integer, intent(in):: Isw
    integer, intent(in):: J
    real(kind=8), intent(in):: Aij(NDEg*NDEg)
    integer, intent(out):: Ir
    !------
    integer:: ndeg2
    integer:: ndeg2l

    Ir = 0
    ndeg2 = NDEg*NDEg
    ndeg2l = NDEg*(NDEg+1)/2
    if ( STAge==30 ) write (6,*) 'warning a matrix was build up '//'but never solved.'
    if ( STAge==10 ) then
      allocate (DIAg(NEQns*ndeg2l),stat=IERror)
      RALoc = RALoc + NEQns*ndeg2l
      allocate (ZLN(LEN_colno*ndeg2),stat=IERror)

      RALoc = RALoc + LEN_colno*ndeg2
      allocate (DSLn(LEN_dsln*ndeg2),stat=IERror)

      if ( IERror/=0 ) stop "Allocation error dsln"

      RALoc = RALoc + LEN_dsln*ndeg2
    endif
    if ( STAge/=20 ) then
      !
      ! for diagonal
      !
      DIAg = 0.
      !
      ! for lower triangle
      !
      ZLN = 0.
      !
      ! for dense window
      !
      DSLn = 0.

      STAge = 20
    endif
    !         Print *,'********Set Stage 20 *********'
    !
    if ( NDEg<=2 ) then
      call ADDR0(Isw,I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,ndeg2,ndeg2l,Ir)
    elseif ( NDEg==3 ) then
      call ADDR3(I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,Ir)
    else
      call ADDRX(I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,NDEg,ndeg2l,Ir)
    endif
  end subroutine STAIJ1

  !======================================================================!
  !> @brief SUM
  !======================================================================!
  subroutine sum(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    implicit none
    !------
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Par(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Diag(*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Zln(*)
    !------
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    real(kind=8):: piv
    real(kind=8):: s
    real(kind=8):: t
    real(kind=8):: zz

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
    integer, intent(in):: Ic
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Zln(*)
    !------
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: ks
    real(kind=8):: s
    real(kind=8):: t
    real(kind=8):: zz

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
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(*)
    real(kind=8), intent(out):: Dsln(*)
    real(kind=8), intent(out):: Temp(*)
    real(kind=8), intent(out):: Zln(*)
    !------
    integer:: ic
    integer:: j
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: ks
    real(kind=8):: s

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
    integer, intent(in):: N
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(*)
    real(kind=8), intent(out):: Dsln(*)
    real(kind=8), intent(out):: Temp(*)
    !------
    integer:: i
    integer:: j
    integer:: joc

    if ( N>0 ) then
      Indx(1) = 0
      joc = 1
      Diag(1) = 1.0D0/Diag(1)
      do i = 2, N
        Indx(i) = joc
        do j = 1, i - 1
          Dsln(joc) = Dsln(joc) - DDOT(Dsln(Indx(i)),Dsln(Indx(j)),j-1)
          joc = joc + 1
        enddo
        call VPROD(Dsln(Indx(i)),Diag,Temp,i-1)
        Diag(i) = Diag(i) - DDOT(Temp,Dsln(Indx(i)),i-1)
        call VCOPY(Temp,Dsln(Indx(i)),i-1)
        Diag(i) = 1.0D0/Diag(i)
      enddo
    endif
  end subroutine SUM3

  !======================================================================!
  !> @brief SXUM
  !======================================================================!
  subroutine SXUM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx,Ndeg,Ndegl,Zz,T)
    implicit none
    !------
    integer, intent(in):: Ic
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(in):: Par(*)
    integer, intent(out):: Indx(*)
    integer, intent(out):: Nch(*)
    real(kind=8), intent(out):: Diag(Ndegl,*)
    real(kind=8), intent(out):: T(Ndegl)
    real(kind=8), intent(out):: Temp(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Zln(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Zz(Ndeg,Ndeg)
    !------
    integer:: ir
    integer:: j
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    integer:: m
    integer:: n
    integer:: ndeg22

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
    integer, intent(in):: Ndeg
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: S(Ndeg,Ndeg)
    real(kind=8), intent(out):: Temp(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Zln(Ndeg,Ndeg,*)
    !------
    integer:: Ic
    integer:: j
    integer:: jc
    integer:: jj
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    integer:: m
    integer:: n

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
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    integer, intent(in):: Neqns
    integer, intent(in):: Nstop
    integer, intent(in):: Xlnzr(*)
    integer, intent(in):: Colno(*)
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(Ndegl,*)
    real(kind=8), intent(out):: Dsln(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Temp(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: Zln(Ndeg,Ndeg,*)
    !------
    integer:: ic
    integer:: j
    integer:: j1
    integer:: j2
    integer:: jc
    integer:: jj
    integer:: joc
    integer:: k
    integer:: ke
    integer:: kk
    integer:: ks
    integer:: locd
    integer:: m
    integer:: n

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
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    integer, intent(in):: Nn
    integer, intent(out):: Indx(*)
    real(kind=8), intent(out):: Diag(Ndegl,*)
    real(kind=8), intent(out):: Dsln(Ndeg,Ndeg,*)
    real(kind=8), intent(out):: T(Ndeg,Ndeg)
    real(kind=8), intent(out):: Temp(Ndeg,Ndeg,*)
    !------
    integer:: i
    integer:: ir
    integer:: j
    integer:: joc
    integer:: locd
    integer:: m
    integer:: n

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
  !> @brief V2PROD
  !======================================================================!
  subroutine V2PROD(A,B,C,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(4,N)
    real(kind=8), intent(in):: B(3,N)
    real(kind=8), intent(out):: C(4,N)
    !------
    integer:: i

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
    integer, intent(in):: N
    real(kind=8), intent(in):: Diag(6,N)
    real(kind=8), intent(in):: Zln(9,N)
    real(kind=8), intent(out):: Zz(9,N)
    !------
    integer:: i

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
  !> @brief VLCPY
  !======================================================================!
  subroutine VLCPY(A,B,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: B(N,N)
    real(kind=8), intent(out):: A(N,N)
    !------
    integer:: i

    do i = 1, N
      A(1:N,i) = B(i,1:N)
    enddo
  end subroutine VLCPY

  !======================================================================!
  !> @brief VCOPY
  !======================================================================!
  subroutine VCOPY(A,C,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(N)
    real(kind=8), intent(out):: C(N)
    !------
    C = A
  end subroutine VCOPY

  !======================================================================!
  !> @brief VPROD
  !======================================================================!
  subroutine VPROD(A,B,C,N)
    implicit none
    !------
    integer, intent(in):: N
    real(kind=8), intent(in):: A(N)
    real(kind=8), intent(in):: B(N)
    real(kind=8), intent(out):: C(N)
    !------
    C(1:N) = A(1:N)*B(1:N)
  end subroutine VPROD

  !======================================================================!
  !> @brief VXPROD
  !======================================================================!
  subroutine VXPROD(Ndeg,Ndegl,Zln,Diag,Zz,N)
    implicit none
    !------
    integer, intent(in):: Ndeg
    integer, intent(in):: Ndegl
    real(kind=8), intent(in):: Diag(Ndegl,N)
    real(kind=8), intent(out):: Zln(Ndeg*Ndeg,N)
    real(kind=8), intent(out):: Zz(Ndeg*Ndeg,N)
    !------
    integer:: i
    integer:: N

    do i = 1, N
      call INVXX(Zz(1,i),Zln(1,i),Diag(1,i),Ndeg)
    enddo
  end subroutine VXPROD
end module HECMW_SOLVER_DIRECT
