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
!> \brief This module contains subroutines of fortran dummy interface to Metis 4.0.3 

module m_Metis403API
  implicit none
!#ifdef USE_64BIT
!  integer,parameter ::  idx_t = c_long_long
!  integer,parameter ::  real_t= c_double
!#else
  integer,parameter ::  idx_t = 4
  integer,parameter ::  real_t= 8
!#endif
!< Return codes
  integer,parameter ::  METIS_OK              = 1   !< Returned normally */
  integer,parameter ::  METIS_ERROR_INPUT     = -2  !< Returned due to erroneous inputs and/or options */
  integer,parameter ::  METIS_ERROR_MEMORY    = -3  !< Returned due to insufficient memory */
  integer,parameter ::  METIS_ERROR           = -4  !< Some other errors */
!
  integer,parameter ::  METIS_ERROR_SETUP     = -11 !< Returned due to setup data contradiction for METIS5
  integer,parameter ::  METIS_ERROR_SUPERVERTEX = -12 !< Returned due to error inputs of super vertices for extension functions
!
!< Metis 4.0.3 partition method ID
  integer,parameter ::  METIS_TYPE_PartGraphRecursive   = 1
  integer,parameter ::  METIS_TYPE_mCPartGraphRecursive = 2
  integer,parameter ::  METIS_TYPE_WPartGraphRecursive  = 3
  
  integer,parameter ::  METIS_TYPE_PartGraphKway        = 1
  integer,parameter ::  METIS_TYPE_PartGraphVKway       = 2
  integer,parameter ::  METIS_TYPE_mCPartGraphKway      = 3
  integer,parameter ::  METIS_TYPE_WPartGraphKway       = 4
  integer,parameter ::  METIS_TYPE_WPartGraphVKway      = 5
  
!<  Extension functions of Metis 5.0.2
!!  & Data structures for fortran codes
  type METIS4
    integer(idx_t)                ::  nvtxs   = 0
    integer(idx_t)                ::  ncon    = 0
    integer(idx_t)                ::  nparts  = 0
    integer(idx_t),pointer        ::  xadj(:)=>NULL(),adjncy(:)=>NULL()
    integer(idx_t),pointer        ::  vwgt(:)=>NULL(),vsize(:)=>NULL(),adjwgt(:)=>NULL()
    real(real_t),pointer          ::  tpwgts(:)=>NULL(),ubvec(:)=>NULL()
!    integer(idx_t),allocatable    ::  vwgt(:),vsize(:),adjwgt(:)
!    real(real_t),allocatable      ::  tpwgts(:),ubvec(:)
    integer(idx_t),pointer        ::  options(:)=>NULL()
    integer(idx_t)                ::  objval = 0
    integer(idx_t),pointer        ::  part(:)=>NULL()
    integer(idx_t)                ::  wgtflag = 0
    integer(idx_t)                ::  optype  = 1
!
!<  Mesh
!    integer(idx_t)                ::  ne  = 0,nn  = 0
!    integer(idx_t),pointer        ::  eptr(:)=>NULL(),eind(:)=>NULL()
!    integer(idx_t)                ::  ncommon = 1
!    integer(idx_t),pointer        ::  epart(:)=>NULL(),npart(:)=>NULL()
!
!<  Reordering
    integer(idx_t),pointer        ::  perm(:)=>NULL(),iperm(:)=>NULL()
!
!<  Super Vertex
    integer(idx_t)                ::  nsvtxs  = 0
    integer(idx_t),pointer        ::  svptr(:)=>NULL()
    integer(idx_t),pointer        ::  svind(:)=>NULL()
    integer(idx_t),pointer        ::  idxOld2New(:)=>NULL()
  end type

  type,private  ::  CONNECT
    integer(idx_t)            ::  num = 0
    integer(idx_t),pointer    ::  idx(:)=>NULL()
  end type
  
  interface freeConnect
    module procedure    freeConnect1
    module procedure    freeConnectM
  end interface
  private ::  freeConnect,freeConnect1,freeConnectM,  &
              MetisF_Init,MetisF_Free

contains
!<  Private function
subroutine freeConnect1(con)
  type(CONNECT),intent(inout) ::  con
!
  integer ::  istat
    if(associated(con%idx)) deallocate(con%idx,stat=istat)
end subroutine  freeConnect1

subroutine freeConnectM(con)
  type(CONNECT),intent(inout) ::  con(:)
!
  integer ::  i,istat
    if(size(con) == 0) return
    do i=1,size(con)
      call freeConnect1(con(i))
    enddo
end subroutine freeConnectM

subroutine MetisF_Init(mts)
  type(METIS4),intent(inout)  ::  mts
!
    mts%nvtxs   = 0
    mts%ncon    = 1
    mts%nparts  = 1
    mts%wgtflag = 0
!
!    mts%ne      = 0
!    mts%nn      = 0
!    mts%ncommon = 1
!
    mts%nsvtxs  = 0
    call MetisF_Free(mts)
end subroutine MetisF_Init

subroutine MetisF_Free(mts)
  type(METIS4),intent(inout)  ::  mts
!
  integer(idx_t)  ::  istat
    if(associated(mts%xadj))    deallocate(mts%xadj,stat=istat)
    if(associated(mts%adjncy))  deallocate(mts%adjncy,stat=istat)
    if(associated(mts%vwgt))    deallocate(mts%vwgt,stat=istat)
    if(associated(mts%vsize))   deallocate(mts%vsize,stat=istat)
    if(associated(mts%adjwgt))  deallocate(mts%adjwgt,stat=istat)
    if(associated(mts%tpwgts))  deallocate(mts%tpwgts,stat=istat)
    if(associated(mts%ubvec))   deallocate(mts%ubvec,stat=istat)
    
    if(associated(mts%options)) deallocate(mts%options,stat=istat)
    if(associated(mts%part))    deallocate(mts%part,stat=istat)
!
    if(associated(mts%perm))  deallocate(mts%perm,stat=istat)
    if(associated(mts%iperm)) deallocate(mts%iperm,stat=istat)
!
    if(associated(mts%svptr))       deallocate(mts%svptr,stat=istat)
    if(associated(mts%svind))       deallocate(mts%svind,stat=istat)
    if(associated(mts%idxOld2New))  deallocate(mts%idxOld2New,stat=istat)    
end subroutine MetisF_Free

!<  Fortran Interface to Metis 4.0.3
function METIS_PartGraphRecursive                   &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  options,edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  edgecut
  integer(idx_t),pointer      ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_PartGraphRecursive

function METIS_PartGraphKway                        &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  options,edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  edgecut
  integer(idx_t),pointer      ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_PartGraphKway

function METIS_PartGraphVKway                       &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   vsize,  wgtflag,  numflag,  &
               nparts,  options,volume,   part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),vsize(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  volume
  integer(idx_t),pointer            ::  part(:)
!
  volume = 0
  ierr = 0
end function METIS_PartGraphVKway

function METIS_mCPartGraphRecursive                 &
              (nvtxs,   ncon,   xadj,     adjncy,   &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  options,edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_mCPartGraphRecursive

function METIS_mCPartGraphKway                      &
              (nvtxs,   ncon,   xadj,     adjncy,   &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  ubvec,  options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  ubvec(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_mCPartGraphKway

function METIS_WPartGraphRecursive                  &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_WPartGraphRecursive

function METIS_WPartGraphKway                       &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  edgecut = 0
  ierr = 0
end function METIS_WPartGraphKway

function METIS_WPartGraphVKway                      &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  volume,   part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  volume
  integer(idx_t),pointer            ::  part(:)
!
  volume = 0
  ierr = 0
end function METIS_WPartGraphVKway

function METIS_EdgeND                             &
              (nvtxs,   xadj,   adjncy, numflag,  &
               options, perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  ierr = 0
end function METIS_EdgeND

function METIS_NodeND                             &
              (nvtxs,   xadj,   adjncy, numflag,  &
               options, perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  ierr = 0
end function METIS_NodeND

function METIS_NodeWND                          &
              (nvtxs,   xadj,   adjncy, vwgt,   &
               numflag, options,perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  ierr = 0
end function METIS_NodeWND

function METIS_EstimateMemory                      &
              (nvtxs,   xadj,   adjncy, numflag,  &
               optype,  nbytes) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag,optype
  integer(idx_t),intent(out)        ::  nbytes
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
!
  ierr = 0
end function METIS_EstimateMemory

!<  Extend METIS_PartGraphRecursive with super vertex information
subroutine MetisF_PartGraphRecursive_Type(mts)
  type(METIS4),intent(inout)  ::  mts
!
    if(associated(mts%tpwgts)) then
      mts%optype = METIS_TYPE_WPartGraphRecursive
      if(.not.associated(mts%vwgt).and.(.not.associated(mts%adjwgt))) then
        mts%wgtflag = 0
      elseif(.not.associated(mts%vwgt)) then
        mts%wgtflag = 1
      elseif(.not.associated(mts%adjwgt)) then
        mts%wgtflag = 2
      else
        mts%wgtflag = 3
      endif
    elseif(mts%ncon == 1) then
      mts%optype = METIS_TYPE_PartGraphRecursive
      if(.not.associated(mts%vwgt).and.(.not.associated(mts%adjwgt))) then
        mts%wgtflag = 0
      elseif(.not.associated(mts%vwgt)) then
        mts%wgtflag = 1
      elseif(.not.associated(mts%adjwgt)) then
        mts%wgtflag = 2
      else
        mts%wgtflag = 3
      endif
    elseif(mts%ncon > 1) then
      mts%optype = METIS_TYPE_mCPartGraphRecursive
      if(.not.associated(mts%adjwgt)) then
        mts%wgtflag = 0
      else
        mts%wgtflag = 1
      endif
    else
      stop 'Error: Unknown type for PartGraphRecursive'
    endif
end subroutine MetisF_PartGraphRecursive_Type

subroutine MetisF_PartGraphKway_Type(mts)
  type(METIS4),intent(inout)  ::  mts
!
    if(associated(mts%tpwgts)) then
      if(associated(mts%adjwgt)) then
        mts%optype = METIS_TYPE_WPartGraphKway
      elseif(associated(mts%vsize)) then
        mts%optype = METIS_TYPE_WPartGraphVKway
      endif
      if( .not.associated(mts%vwgt).and.      &
         (.not.associated(mts%adjwgt)).and.   &
         (.not.associated(mts%vsize))) then
        mts%wgtflag = 0
      elseif(.not.associated(mts%vwgt).and.   &
         (associated(mts%adjwgt).or.          &
          associated(mts%vsize))) then
        mts%wgtflag = 1
      elseif(associated(mts%vwgt).and.        &
         (.not.associated(mts%adjwgt)).and.   &
         (.not.associated(mts%vsize))) then
        mts%wgtflag = 2
      else
        mts%wgtflag = 3
      endif
    elseif(mts%ncon > 1.and.associated(mts%ubvec)) then
      mts%optype = METIS_TYPE_mCPartGraphKway
      if(.not.associated(mts%adjwgt)) then
        mts%wgtflag = 0
      else
        mts%wgtflag = 1
      endif
    elseif(.not.associated(mts%ubvec)) then
      if(associated(mts%adjwgt)) then
        mts%optype = METIS_TYPE_PartGraphKway
      elseif(associated(mts%vsize)) then
        mts%optype = METIS_TYPE_PartGraphVKway
      endif
      if( .not.associated(mts%vwgt).and.      &
         (.not.associated(mts%adjwgt)).and.   &
         (.not.associated(mts%vsize))) then
        mts%wgtflag = 0
      elseif(.not.associated(mts%vwgt).and.   &
         (associated(mts%adjwgt).or.          &
          associated(mts%vsize))) then
        mts%wgtflag = 1
      elseif(associated(mts%vwgt).and.        &
         (.not.associated(mts%adjwgt)).and.   &
         (.not.associated(mts%vsize))) then
        mts%wgtflag = 2
      else
        mts%wgtflag = 3
      endif
    else
      stop 'Error: Unknown type for PartGraphKway'
    endif
end subroutine MetisF_PartGraphKway_Type

end module m_Metis403API
