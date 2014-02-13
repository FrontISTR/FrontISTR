!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Jun Yin (ASTOM)                                !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines to internal tool mesh structure
!   for parallel contact simulation

module m_Makrose
!  use m_Kind
  use hecmw
  implicit none
! Mesh types
  integer,parameter     ::  MAKROSE_POINT = 1
  integer,parameter     ::  MAKROSE_LINE  = 20
  integer,parameter     ::  MAKROSE_LINE3 = 22
  integer,parameter     ::  MAKROSE_LINE4 = 23
!
  integer,parameter     ::  MAKROSE_TRI3  = 30
  integer,parameter     ::  MAKROSE_TRI6  = 32
  integer,parameter     ::  MAKROSE_TRI9  = 33
  integer,parameter     ::  MAKROSE_TRI10 = 37 
  integer,parameter     ::  MAKROSE_QUAD4 = 40
  integer,parameter     ::  MAKROSE_QUAD8 = 42
  integer,parameter     ::  MAKROSE_QUAD12= 43
  integer,parameter     ::  MAKROSE_QUAD9 = 46
  integer,parameter     ::  MAKROSE_QUAD16= 47
!
  integer,parameter     ::  MAKROSE_PRI6  = 60
  integer,parameter     ::  MAKROSE_PRI15 = 62
  integer,parameter     ::  MAKROSE_TET4  = 70
  integer,parameter     ::  MAKROSE_TET10 = 72
  integer,parameter     ::  MAKROSE_TET16 = 73
  integer,parameter     ::  MAKROSE_HEXA8 = 80
  integer,parameter     ::  MAKROSE_HEXA20= 82
  integer,parameter     ::  MAKROSE_HEXA32= 83
!
  integer,parameter     ::  MAKROSE_PYRD5 = 90
  integer,parameter     ::  MAKROSE_PYRD13= 92
! 
  integer               ::  MAKROSE_NODE_NUM(99)
  data MAKROSE_NODE_NUM / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 2,   &
                          0, 3, 4, 0, 0, 0, 0, 0, 0, 3,   &
                          0, 6, 9, 0, 0, 0,10, 0, 0, 4,   &
                          0, 8,12, 0, 0, 9,16, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 6,   &
                          0,15, 0, 0, 0, 0, 0, 0, 0, 4,   &
                          0,10,16, 0, 0, 0, 0, 0, 0, 8,   &
                          0,20,32, 0, 0, 0, 0, 0, 0, 5,   &
                          0,13, 0, 0, 0, 0, 0, 0, 0/
  integer               ::  MAKROSE_EDGE_NUM(99)
  data MAKROSE_EDGE_NUM / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   &
                          0, 1, 1, 0, 0, 0, 0, 0, 0, 3,   &
                          0, 3, 3, 0, 0, 0, 3, 0, 0, 4,   &
                          0, 4, 4, 0, 0, 4, 4, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 9,   &
                          0, 9, 0, 0, 0, 0, 0, 0, 0, 6,   &
                          0, 6, 6, 0, 0, 0, 0, 0, 0,12,   &
                          0,12,12, 0, 0, 0, 0, 0, 0, 8,   &
                          0, 8, 0, 0, 0, 0, 0, 0, 0/
  integer               ::  MAKROSE_FACE_NUM(99)
  data MAKROSE_FACE_NUM / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   &
                          0, 1, 1, 0, 0, 0, 1, 0, 0, 1,   &
                          0, 1, 1, 0, 0, 1, 1, 0, 0, 0,   &
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 5,   &
                          0, 5, 0, 0, 0, 0, 0, 0, 0, 4,   &
                          0, 4, 4, 0, 0, 0, 0, 0, 0, 6,   &
                          0, 6, 6, 0, 0, 0, 0, 0, 0, 5,   &
                          0, 5, 0, 0, 0, 0, 0, 0, 0/
  integer               ::  MAKROSE_EDGENODE_NUM(99)
  data MAKROSE_EDGENODE_NUM / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 4, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 4, 0, 0, 0, 4, 0, 0, 2,   &
                              0, 3, 4, 0, 0, 3, 4, 0, 0, 0,   &
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 0, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 4, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 4, 0, 0, 0, 0, 0, 0, 2,   &
                              0, 3, 0, 0, 0, 0, 0, 0, 0/
  include "MeshDefinition.inc"
!
  type MAKROSE_STRUCT
    integer(kint)           ::  ndim = 3
    integer(kint)           ::  nn = 0
    integer(kint)           ::  nn_i = 0
    integer(kint),pointer   ::  ngid(:)=>null()   ! Nodal global ID
    real(kreal),pointer     ::  x(:,:)=>null()
    integer(kint)           ::  ne = 0
    integer(kint)           ::  ne_i = 0
    integer(kint),pointer   ::  egid(:)=>null()   ! Element global ID
    integer(kint),pointer   ::  eptr(:)=>null()   ! Elemental starting node pointer
    integer(kint),pointer   ::  eind(:)=>null()   ! Elemental node index
    integer(kint),pointer   ::  etyp(:)=>null()   ! Element type
    integer(kint),pointer   ::  egrp(:)=>null()   ! Element group ID
    integer(kint),pointer   ::  emat(:)=>null()   ! Element material ID
!
    integer(kint),pointer   ::  nptr(:)=>null()   ! Number of elements that share a node
    integer(kint),pointer   ::  nind(:)=>null()   ! Element index sharing a node
  end type
  
!  interface size
!    module procedure sizeNonzero
!  end interface

contains

function sizeNonzero(iarray) result(n)
  integer,intent(in)  ::  iarray(:)
  integer   ::  n
!
  integer ::  i
    do i=1,size(iarray)
      if(iarray(i) == 0) exit
    enddo
    n = i - 1
end function sizeNonzero    

function Mak_GetMeshNodeToNode(mtype,idx) result(nidx)
  integer,intent(in)  ::  mtype
  integer,intent(in)  ::  idx
  integer             ::  nidx(6)
!
    nidx(:) = 0
    select case(mtype)
    case(MAKROSE_POINT)
    case(MAKROSE_LINE)
      nidx(1) = LINE%n_n1(1,idx)
    case(MAKROSE_LINE3)
      nidx(1:2) = LINE%n_n1(1:2,idx)
    case(MAKROSE_LINE4 )
      stop
    case(MAKROSE_TRI3  )
      nidx(1:2) = TRI%n_n1(1:2,idx)
    case(MAKROSE_TRI6  )
      nidx(1:4) = TRI%n_n1(1:4,idx)
    case(MAKROSE_TRI9  )
      stop
    case(MAKROSE_TRI10 )
      stop
    case(MAKROSE_QUAD4 )
      nidx(1:2) = QUAD%n_n1(1:2,idx)
    case(MAKROSE_QUAD8 )
      nidx(1:4) = QUAD%n_n1(1:4,idx)
    case(MAKROSE_QUAD12)
      stop
    case(MAKROSE_QUAD9 )
      stop
    case(MAKROSE_QUAD16)
      stop
    case(MAKROSE_PRI6  )
      nidx(1:3) = PRISM%n_n1(1:3,idx)
    case(MAKROSE_PRI15 )
      nidx(1:6) = PRISM%n_n1(1:6,idx)
    case(MAKROSE_TET4  )
      nidx(1:3) = TETRA%n_n1(1:3,idx)
    case(MAKROSE_TET10 )
      nidx(1:6) = TETRA%n_n1(1:6,idx)
    case(MAKROSE_TET16 )
      stop
    case(MAKROSE_HEXA8 )
      nidx(1:3) = HEXA%n_n1(1:3,idx)
    case(MAKROSE_HEXA20)
      nidx(1:6) = HEXA%n_n1(1:6,idx)
    case(MAKROSE_HEXA32)
      stop
    case(MAKROSE_PYRD5 )
      if(idx <= 4) then
        nidx(1:3) = PYRAMID%n_n1(1:3,idx)
      else
        nidx(1:4) = PYRAMID%n5_n1(1:4)
      endif
    case(MAKROSE_PYRD13)
      stop
    case default
      stop
    end select
end function Mak_GetMeshNodeToNode
    
subroutine Mak_Init(mak,ndim,nn,ne)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
  integer(kint),intent(in),optional   ::  ndim
  integer(kint),intent(in),optional   ::  nn
  integer(kint),intent(in),optional   ::  ne
!
  integer   ::  istat
    mak%ndim  = 3
    mak%nn    = 0;  mak%nn_i  = 0
    if(associated(mak%ngid)) deallocate(mak%ngid,stat=istat)
    if(associated(mak%x)) deallocate(mak%x,stat=istat)
    mak%ne    = 0;  mak%ne_i  = 0
    if(associated(mak%egid)) deallocate(mak%egid,stat=istat)
    if(associated(mak%eptr)) deallocate(mak%eptr,stat=istat)
    if(associated(mak%eind)) deallocate(mak%eind,stat=istat)
    if(associated(mak%etyp)) deallocate(mak%etyp,stat=istat)
    if(associated(mak%egrp)) deallocate(mak%egrp,stat=istat)
    if(associated(mak%emat)) deallocate(mak%emat,stat=istat)
!
    if(associated(mak%nptr)) deallocate(mak%nptr,stat=istat)
    if(associated(mak%nind)) deallocate(mak%nind,stat=istat)
!
    if(present(ndim)) mak%ndim = ndim
    if(present(nn)) then
      mak%nn = nn
      allocate(mak%ngid(nn),stat=istat)
      allocate(mak%x(mak%ndim,mak%nn),stat=istat)
    endif
    if(present(ne)) then
      mak%ne = ne
      allocate(mak%egid(ne),stat=istat)
      allocate(mak%etyp(ne),stat=istat)
      allocate(mak%egrp(ne),stat=istat)
      allocate(mak%emat(ne),stat=istat)
    endif     
end subroutine Mak_Init

subroutine Mak_ReadFile_tmp(mak,filename)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
  character(len=*),intent(in)         ::  filename
!
  character(len=256)  ::  head,filename_tmp
  logical   ::  log
  integer   ::  i,unit=13,n,istat,maxGID=0,idummy(20)
  integer,pointer ::  nlid(:)=>null(),eind(:)=>null()
!
    n = len_trim(filename)
    if(n > 4) then
      if(filename(n-3:n) == '.efp'.or.  &
         filename(n-3:n) == '.efe') then
        head = filename(1:n-4)
        n = n - 4
      else
        head = filename
      endif
    else
      head = filename
    endif
!   Read nodal coordinates
    filename_tmp = head(1:n)//'.efp'
    inquire(file=filename_tmp,exist=log)
    if(.not.log) then
      print *,'file ',filename_tmp(1:len_trim(filename_tmp)),' does not exist!'
      stop
    endif
    open(unit,file=filename_tmp,status='unknown')
    print *,'Enter Node Number:'
    read(*,*)mak%nn
    allocate(mak%x(3,mak%nn),stat=istat)
    allocate(mak%ngid(mak%nn),stat=istat)
    do i=1,mak%nn
      read(unit,*)mak%ngid(i),mak%x(1:3,i)
      maxGID = max(mak%ngid(i),maxGID)
    enddo
    close(unit)
    allocate(nlid(maxGID),stat=istat)
    nlid(:) = 0
    do i=1,mak%nn
      nlid(mak%ngid(i)) = i
    enddo
!   Read element data
    filename_tmp = head(1:n)//'.efe'
    inquire(file=filename_tmp,exist=log)
    if(.not.log) then
      print *,'file ',filename_tmp(1:len_trim(filename_tmp)),' does not exist!'
      stop
    endif
    open(unit,file=filename_tmp,status='unknown')
    print *,'Enter Element Number:'
    read(*,*)mak%ne
    allocate(mak%eptr(mak%ne+1),stat=istat)
    allocate(mak%egid(mak%ne),stat=istat)
    allocate(mak%etyp(mak%ne),stat=istat)
    allocate(mak%emat(mak%ne),stat=istat)
    allocate(mak%egrp(mak%ne),stat=istat)
    allocate(eind(mak%ne*20),stat=istat)
    mak%eptr(1) = 1
    do i=1,mak%ne
      read(unit,*)mak%egid(i),mak%etyp(i),mak%emat(i),mak%egrp(i)
      backspace(unit)
      mak%eptr(i+1) = mak%eptr(i) + MAKROSE_NODE_NUM(mak%etyp(i))
      if(mak%eptr(i+1)-1 > size(eind)) then
        stop
      endif
      read(unit,*)idummy(1:4),eind(mak%eptr(i):mak%eptr(i+1)-1)
    enddo
    close(unit)
    allocate(mak%eind(mak%eptr(mak%ne+1)-1),stat=istat)
    do i=1,mak%eptr(mak%ne+1)-1
      mak%eind(i) = nlid(eind(i))
    enddo
    deallocate(nlid,stat=istat)
    deallocate(eind,stat=istat)
end subroutine Mak_ReadFile_tmp

subroutine Mak_ReadFile(mak,filename)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
  character(len=*),intent(in)         ::  filename
!
  character(len=256)  ::  head,filename_tmp
  logical   ::  log
  integer   ::  unit=14,n
!
    n = len_trim(filename)
    if(n > 4) then
      if(filename(n-3:n) == '.efp'.or.  &
         filename(n-3:n) == '.efe') then
        head = filename(1:n-4)
        n = n - 4
      else
        head = filename
      endif
    else
      head = filename
    endif
!   Read nodal coordinates
    filename_tmp = head(1:n)//'.efp'
    inquire(file=filename_tmp,exist=log)
    if(.not.log) then
      print *,'file ',filename_tmp(1:len_trim(filename_tmp)),' does not exist!'
      stop
    endif
    open(unit,file=filename_tmp,status='unknown')
    close(unit)
end subroutine Mak_ReadFile

subroutine Mak_WriteFile(mak,filename,u,skipElmtID)
  type(MAKROSE_STRUCT),intent(in)   ::  mak
  character(len=*),intent(in)       ::  filename
  real(kreal),intent(in),optional   ::  u(:)
  integer(kint),pointer,optional    ::  skipElmtID(:)
!
  character(len=256)  ::  head,filename_tmp
  logical   ::  log
  integer   ::  i,unit=14,n,istat
  integer(kint),allocatable   ::  help(:)
!
    n = len_trim(filename)
    if(n > 4) then
      if(filename(n-3:n) == '.efp'.or.  &
         filename(n-3:n) == '.efe') then
        head = filename(1:n-4)
        n = n - 4
      else
        head = filename
      endif
    else
      head = filename
    endif
    filename_tmp = head(1:n)//'.efp'
    open(unit,file=filename_tmp,status='unknown')
    do i=1,mak%nn
      if(present(u)) then
        write(unit,'(I10,3ES16.8)')i,mak%x(1:3,i)+u(3*(i-1)+1:3*i)
      else
        write(unit,'(I10,3ES16.8)')i,mak%x(1:3,i)
      endif
    enddo
    close(unit)

    if(present(skipElmtID)) then
      if(associated(skipElmtID)) then
        allocate(help(mak%ne),stat=istat)
        help(:) = 0
        help(skipElmtID(:)) = 1
      endif
    endif
!   Read element data
    filename_tmp = head(1:n)//'.efe'
    open(unit,file=filename_tmp,status='unknown')
    do i=1,mak%ne
      if(present(skipElmtID)) then
        if(help(i) == 1) cycle
      endif
      write(unit,'(30I9)')i,mak%etyp(i),mak%emat(i),mak%egrp(i),mak%eind(mak%eptr(i):mak%eptr(i+1)-1)
    enddo
    close(unit)
end subroutine Mak_WriteFile

subroutine Mak_WriteFile_box(mak,xc,dx,ds_x,filename)
  type(MAKROSE_STRUCT),intent(in)   ::  mak
  real(kreal),intent(in)            ::  xc(3),dx(3),ds_x(3)
  character(len=*),intent(in)       ::  filename
!
  character(len=256)  ::  head,filename_tmp
  logical   ::  log
  integer   ::  i,j,k,l,unit=14,n,nn
!
    n = len_trim(filename)
    nn = n
    if(n > 4) then
      if(filename(n-3:n) == '.efp'.or.  &
         filename(n-3:n) == '.efe') then
        head = filename(1:n-4)
        n = n - 4
      else
        head = filename
      endif
    else
      head = filename
    endif
    filename_tmp = head(1:n)//'.efp'
    open(unit,file=filename_tmp,status='unknown')
    do i=1,mak%nn
      write(unit,'(I10,3ES16.8)')i,mak%x(1:3,i)
    enddo
    n = 0
    do l=1,2
    do k=1,2
      do j=1,2
        do i=1,2
          n = n + 1
          if(l == 1) then
            write(unit,'(I10,3ES16.8)')mak%nn+n,xc(1)+(-1.0D0)**i*dx(1),  &
                                                xc(2)+(-1.0D0)**j*dx(2),  &
                                                xc(3)+(-1.0D0)**k*dx(3)
          else
            write(unit,'(I10,3ES16.8)')mak%nn+n,xc(1)+(-1.0D0)**i*ds_x(1),  &
                                                xc(2)+(-1.0D0)**j*ds_x(2),  &
                                                xc(3)+(-1.0D0)**k*ds_x(3)
          endif
        enddo
      enddo
    enddo
    enddo
    
    close(unit)
    
!   Read element data
    filename_tmp = head(1:nn)//'.efe'
    open(unit,file=filename_tmp,status='unknown')
    do i=1,mak%ne
      write(unit,'(30I9)')i,mak%etyp(i),mak%emat(i),mak%egrp(i),mak%eind(mak%eptr(i):mak%eptr(i+1)-1)
    enddo
    n = 0
    do i=1,2
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+1,mak%nn+(i-1)*8+1
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+3,mak%nn+(i-1)*8+4
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+5,mak%nn+(i-1)*8+6
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+7,mak%nn+(i-1)*8+8
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+1,mak%nn+(i-1)*8+3
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+2,mak%nn+(i-1)*8+4
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+5,mak%nn+(i-1)*8+7
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+6,mak%nn+(i-1)*8+8
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+1,mak%nn+(i-1)*8+5
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+2,mak%nn+(i-1)*8+6
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+3,mak%nn+(i-1)*8+7
      n = n + 1
      write(unit,'(30I9)')mak%ne+n,MAKROSE_LINE,1,20+i,mak%nn+(i-1)*8+4,mak%nn+(i-1)*8+8
    enddo
    close(unit)
end subroutine Mak_WriteFile_box

subroutine Mak_GetElementOnNode(mak)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
!
  integer(kint)   ::  i,j,istat
  integer(kint),allocatable   ::  en(:)
!
    allocate(en(mak%nn),stat=istat)
    en(:) = 0
    do i=1,mak%ne
      do j=mak%eptr(i),mak%eptr(i+1)-1
        en(mak%eind(j)) = en(mak%eind(j)) + 1
      enddo
    enddo
    if(associated(mak%nptr)) deallocate(mak%nptr,stat=istat)
    allocate(mak%nptr(mak%nn+1),stat=istat)
    mak%nptr(1) = 1
    do i=1,mak%nn
      mak%nptr(i+1) = mak%nptr(i) + en(i)
    enddo
    if(associated(mak%nind)) deallocate(mak%nind,stat=istat)
    allocate(mak%nind(mak%nptr(mak%nn+1)-1),stat=istat)
    en(:) = 0
    do i=1,mak%ne
      do j=mak%eptr(i),mak%eptr(i+1)-1
        en(mak%eind(j)) = en(mak%eind(j)) + 1
        mak%nind(mak%nptr(mak%eind(j)) - 1 + en(mak%eind(j))) = i
      enddo
    enddo
    deallocate(en,stat=istat)
end subroutine Mak_GetElementOnNode

subroutine Mak_MeshToNodal(mak,xadj,adjncy)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
  integer(kint),pointer               ::  xadj(:)
  integer(kint),pointer               ::  adjncy(:)
!
  integer(kint)   ::  i,j,k,l,n,istat,nid(6),elmtID,nodeID,nn,nnn
  integer(kint),allocatable   ::  help(:)
!
    if(.not.(associated(mak%nptr).and.associated(mak%nind))) then
      call Mak_GetElementOnNode(mak)
    endif
    allocate(help(mak%nn),stat=istat)
    nn = 0
    do i=1,mak%nn
      help(:) = 0
      nnn = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        n = 0
        do k=mak%eptr(elmtID),mak%eptr(elmtID+1)-1
          n = n + 1
          nodeID = mak%eind(k)
          if(nodeID == i) then
            nid(:) = Mak_GetMeshNodeToNode(mak%etyp(elmtID),n)
            do l=1,sizeNonzero(nid)
              if(help(mak%eind(mak%eptr(elmtID)-1+nid(l))) == 0) then
                nnn = nnn + 1
                help(mak%eind(mak%eptr(elmtID)-1+nid(l))) = 1
              endif
            enddo
            exit
          endif
        enddo
      enddo
      nn = nn + nnn
    enddo
!
    allocate(xadj(mak%nn+1),stat=istat)
    xadj(1) = 1
    allocate(adjncy(nn),stat=istat)
!
    nn = 0
    do i=1,mak%nn
      help(:) = 0
      nnn = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        n = 0
        do k=mak%eptr(elmtID),mak%eptr(elmtID+1)-1
          n = n + 1
          nodeID = mak%eind(k)
          if(nodeID == i) then
            nid(:) = Mak_GetMeshNodeToNode(mak%etyp(elmtID),n)
            do l=1,sizeNonzero(nid)
              if(help(mak%eind(mak%eptr(elmtID)-1+nid(l))) == 0) then
                nnn = nnn + 1
                help(mak%eind(mak%eptr(elmtID)-1+nid(l))) = 1
                adjncy(nn+nnn) = mak%eind(mak%eptr(elmtID)-1+nid(l))
              endif
            enddo
            exit
          endif
        enddo
      enddo
      nn = nn + nnn
      xadj(i+1) = xadj(i) + nnn
    enddo
    deallocate(help,stat=istat)
end subroutine Mak_MeshToNodal

subroutine Mak_MeshToDual(mak,nCommon,xadj,adjncy)
  type(MAKROSE_STRUCT),intent(inout)  ::  mak
  integer(kint),intent(in)            ::  nCommon   ! Number of common nodes between elements
  integer(kint),pointer               ::  xadj(:)
  integer(kint),pointer               ::  adjncy(:)
!
  integer(kint)   ::  i,j,k,l,n,istat,elmtID,nodeID,ne,nne
  integer(kint),allocatable   ::  help(:),counter(:)
!
    
    if(.not.(associated(mak%nptr).and.associated(mak%nind))) then
      call Mak_GetElementOnNode(mak)
    endif
    allocate(help(mak%ne),stat=istat)
    allocate(counter(mak%ne),stat=istat)
    ne = 0
    do i=1,mak%ne
      counter(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=mak%nptr(nodeID),mak%nptr(nodeID+1)-1
          elmtID = mak%nind(k)
          if(elmtID == i) cycle
          counter(elmtID) = counter(elmtID) + 1
        enddo
      enddo
      help(:) = 0
      nne = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=mak%nptr(nodeID),mak%nptr(nodeID+1)-1
          elmtID = mak%nind(k)
          if(elmtID == i) cycle
          if(help(elmtID) /= 0) cycle
          help(elmtID) = 1
          if(counter(elmtID) >= nCommon) then
            nne = nne + 1
          endif
        enddo
      enddo
      ne = ne + nne
    enddo
!
    allocate(xadj(mak%ne+1),stat=istat)
    xadj(1) = 1
    allocate(adjncy(ne),stat=istat)
!
    ne = 0
    do i=1,mak%ne
      counter(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=mak%nptr(nodeID),mak%nptr(nodeID+1)-1
          elmtID = mak%nind(k)
          if(elmtID == i) cycle
          counter(elmtID) = counter(elmtID) + 1
        enddo
      enddo
      help(:) = 0
      nne = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=mak%nptr(nodeID),mak%nptr(nodeID+1)-1
          elmtID = mak%nind(k)
          if(elmtID == i) cycle
          if(help(elmtID) /= 0) cycle
          help(elmtID) = 1        
          if(counter(elmtID) >= nCommon) then
            nne = nne + 1
            adjncy(ne+nne) = elmtID
          endif
        enddo
      enddo
      ne = ne + nne
      xadj(i+1) = xadj(i) + nne
    enddo
    deallocate(help,counter,stat=istat)
end subroutine Mak_MeshToDual        

function Mak_GetElementCommonNodeNumber(mak,e1,e2) result(n)
  type(MAKROSE_STRUCT),intent(in) ::  mak
  integer(kint),intent(in)        ::  e1,e2   ! two element ID
  integer(kint)   ::  n
!
  integer(kint)   ::  i,j,istat,nodeID
  integer(kint),allocatable   ::  help(:)
!
    n = 0
    if(e1 == e2) then
      n = mak%eptr(e1+1) - mak%eptr(e1)
      return
    endif
    allocate(help(mak%nn),stat=istat)
    help(:) = 0
    do i=mak%eptr(e1),mak%eptr(e1+1)-1
      nodeID = mak%eind(i)
      if(help(nodeID) /= 0) cycle
      help(nodeID) = 1
      do j=mak%nptr(nodeID),mak%nptr(nodeID+1)-1
        if(mak%nind(j) == e2) then
          n = n + 1
          exit
        endif
      enddo
    enddo
    deallocate(help,stat=istat)
end function Mak_GetElementCommonNodeNumber
        
subroutine Mak_GetElementInBox(mak,boxCenter,dxyz,n,svind,flagNode)
  type(MAKROSE_STRUCT),intent(in) ::  mak
  real(kreal),intent(in)          ::  boxCenter(:)
  real(kreal),intent(in)          ::  dxyz(:)
  integer(kint),intent(out)       ::  n
  integer(kint),pointer           ::  svind(:)
  logical,intent(in),optional     ::  flagNode
!
  integer(kint) ::  i,j,k,istat,ndim,nodeID,count
  real(kreal),allocatable   ::  x0(:),dx(:)
  integer(kint),allocatable ::  minID(:),maxID(:),pcID(:),idx(:),help(:)
!
    n = 0
    ndim = min(mak%ndim,size(boxCenter),size(dxyz))
    if(ndim <= 0) stop
    allocate(x0(ndim),stat=istat)
    allocate(dx(ndim),stat=istat)
    x0(1:ndim) = boxCenter(1:ndim) - dxyz(1:ndim)
    dx(1:ndim) = 2.0D0*dxyz(1:ndim)
    allocate(minID(ndim),maxID(ndim),stat=istat)
    allocate(pcID(ndim),stat=istat)
    allocate(idx(mak%ne),stat=istat)
    next_elmt : do i=1,mak%ne
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=1,ndim
          pcID(k) = (mak%x(k,nodeID)-x0(k))/dx(k) + 1
          minID(k) = min(minID(k),pcID(k))
          maxID(k) = max(maxID(k),pcID(k))
        enddo
        if(j == mak%eptr(i)) then
          minID(:) = pcID(:)
          maxID(:) = pcID(:)
        endif
      enddo
      do j=1,ndim
        if(minID(j) > 1.or.maxID(j) < 1) cycle next_elmt
      enddo
      if(minID(1) < -1.or.maxID(1) > 2) then
!       exclude elements on dummy block boundary
        continue
        cycle next_elmt
      endif
      n = n + 1
      idx(n) = i
    enddo next_elmt
    if(associated(svind)) deallocate(svind,stat=istat)
    if(.not.present(flagNode)) then
      allocate(svind(n),stat=istat)
      svind(1:n) = idx(1:n)
    else
      allocate(help(mak%nn),stat=istat)
      help(:) = 0
      count = 0
      do i=1,n
        do j=mak%eptr(idx(i)),mak%eptr(idx(i)+1)-1
          if(help(mak%eind(j)) == 0) then
            count = count + 1
            help(mak%eind(j)) =1
          endif
        enddo
      enddo
      allocate(svind(count),stat=istat)
      count = 0
      do i=1,mak%nn
        if(help(i) == 0) cycle
        count = count + 1
        svind(count) = i
      enddo
      deallocate(help,stat=istat)
      n = count
    endif   
    deallocate(idx,x0,dx,minID,maxID,pcID,stat=istat)
end subroutine Mak_GetElementInBox

subroutine Mak_GetElementInBoxEndPoint(mak,boxCenter,dxyz,n,svind,flagNode)
  type(MAKROSE_STRUCT),intent(in) ::  mak
  real(kreal),intent(in)          ::  boxCenter(:)
  real(kreal),intent(in)          ::  dxyz(:)
  integer(kint),intent(out)       ::  n
  integer(kint),pointer           ::  svind(:)
  logical,intent(in),optional     ::  flagNode
!
  integer(kint) ::  i,j,k,istat,ndim,nodeID,count
  real(kreal),allocatable   ::  x0(:),dx(:)
  integer(kint),allocatable ::  minID(:),maxID(:),pcID(:),idx(:),help(:)
!
    n = 0
    ndim = min(mak%ndim,size(boxCenter),size(dxyz))
    if(ndim <= 0) stop
    allocate(x0(ndim),stat=istat)
    allocate(dx(ndim),stat=istat)
    x0(1:ndim) = boxCenter(1:ndim)
    dx(1:ndim) = dxyz(1:ndim)
    allocate(minID(ndim),maxID(ndim),stat=istat)
    allocate(pcID(ndim),stat=istat)
    allocate(idx(mak%ne),stat=istat)
    next_elmt : do i=1,mak%ne
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        do k=1,ndim
          pcID(k) = (mak%x(k,nodeID)-x0(k))/dx(k) + 1
          minID(k) = min(minID(k),pcID(k))
          maxID(k) = max(maxID(k),pcID(k))
        enddo
        if(j == mak%eptr(i)) then
          minID(:) = pcID(:)
          maxID(:) = pcID(:)
        endif
      enddo
      do j=1,ndim
        if(minID(j) > 1.or.maxID(j) < 1) cycle next_elmt
      enddo
      if(minID(1) < -1.or.maxID(1) > 2) then
!       exclude elements on dummy block boundary
        continue
        cycle next_elmt
      endif
      n = n + 1
      idx(n) = i
    enddo next_elmt
    if(associated(svind)) deallocate(svind,stat=istat)
    if(.not.present(flagNode)) then
      allocate(svind(n),stat=istat)
      svind(1:n) = idx(1:n)
    else
      allocate(help(mak%nn),stat=istat)
      help(:) = 0
      count = 0
      do i=1,n
        do j=mak%eptr(idx(i)),mak%eptr(idx(i)+1)-1
          if(help(mak%eind(j)) == 0) then
            count = count + 1
            help(mak%eind(j)) =1
          endif
        enddo
      enddo
      allocate(svind(count),stat=istat)
      count = 0
      do i=1,mak%nn
        if(help(i) == 0) cycle
        count = count + 1
        svind(count) = i
      enddo
      deallocate(help,stat=istat)
      n = count
    endif   
    deallocate(idx,x0,dx,minID,maxID,pcID,stat=istat)
end subroutine Mak_GetElementInBoxEndPoint


end module m_Makrose
