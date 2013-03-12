!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
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
!> \brief This module contains subroutines of fortran interface to Metis 4.0.3

module m_Metis403API
  use, intrinsic  ::  iso_c_binding   !, only : C_int, C_double, C_float, 
  implicit none
!  integer,parameter ::  kint  = kind(10000000000000)
!  integer,parameter	::  kreal = kind(1.1111111111D0)
#ifdef USE_64BIT
  integer,parameter ::  idx_t = c_long_long
  integer,parameter ::  real_t= c_double
#else
  integer,parameter ::  idx_t = c_int
  integer,parameter ::  real_t= c_float
#endif
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

!< Interface definition for calling c functions in Metis 4.0.3
  interface
!<  Graph partitioning routines
    function METIS_PartGraphRecursive_f                 &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  options,edgecut,  part)     &
      result(ierr) bind(C,name="METIS_PartGraphRecursive")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_PartGraphRecursive_f

    function METIS_PartGraphKway_f                      &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  options,edgecut,  part)     &
      result(ierr) bind(C,name="METIS_PartGraphKway")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_PartGraphKway_f

    function METIS_PartGraphVKway_f                     &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   vsize,  wgtflag,  numflag,  &
                   nparts,  options,volume,  part)     &
      result(ierr) bind(C,name="METIS_PartGraphVKway")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,vsize
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  volume
      type(c_ptr),value           ::  part
    end function METIS_PartGraphVKway_f

    function METIS_mCPartGraphRecursive_f               &
                  (nvtxs,   ncon,   xadj,     adjncy,   &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  options,edgecut,  part)     &
      result(ierr) bind(C,name="METIS_mCPartGraphRecursive")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,ncon,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_mCPartGraphRecursive_f

    function METIS_mCPartGraphKway_f                    &
                  (nvtxs,   ncon,   xadj,     adjncy,   &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  ubvec,  options,  edgecut,  part)     &
      result(ierr) bind(C,name="METIS_mCPartGraphKway")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,ncon,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  ubvec
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_mCPartGraphKway_f

    function METIS_WPartGraphRecursive_f                &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  tpwgts, options,edgecut,  part)     &
      result(ierr) bind(C,name="METIS_WPartGraphRecursive")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  tpwgts
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_WPartGraphRecursive_f

    function METIS_WPartGraphKway_f                     &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  tpwgts, options,  edgecut,  part)     &
      result(ierr) bind(C,name="METIS_WPartGraphKway")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  tpwgts
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  edgecut
      type(c_ptr),value           ::  part
    end function METIS_WPartGraphKway_f
    
    function METIS_WPartGraphVKway_f                    &
                  (nvtxs,   xadj,   adjncy,             &
                    vwgt,   adjwgt, wgtflag,  numflag,  &
                   nparts,  tpwgts, options,  volume,   part)     &
      result(ierr) bind(C,name="METIS_WPartGraphVKway")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,nparts
      integer(idx_t),intent(in)   ::  wgtflag,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt,adjwgt
      type(c_ptr),value           ::  tpwgts
      type(c_ptr),value           ::  options
      integer(idx_t),intent(out)  ::  volume
      type(c_ptr),value           ::  part
    end function METIS_WPartGraphVKway_f

!<  Mesh partitioning routines
!!  Because Metis 4.0.3 does not support mixed element type, 
!!  its APIs are not so much useful and the API interfaces are not provided.
!!  Users should provide, by themselves, graph data of mixed element-type mesh and use
!!  graph partition APIs directly

!<  Sparse matrix reordering rountines
    function METIS_EdgeND_f                           &
                  (nvtxs,   xadj,   adjncy, numflag,  &
                   options, perm,   iperm)            &
      result(ierr) bind(C,name="METIS_EdgeND")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  options
      type(c_ptr),value           ::  perm
      type(c_ptr),value           ::  iperm
    end function METIS_EdgeND_f

    function METIS_NodeND_f                           &
                  (nvtxs,   xadj,   adjncy, numflag,  &
                   options, perm,   iperm)            &
      result(ierr) bind(C,name="METIS_NodeND")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  options
      type(c_ptr),value           ::  perm
      type(c_ptr),value           ::  iperm
    end function METIS_NodeND_f

    function METIS_NodeWND_f                        &
                  (nvtxs,   xadj,   adjncy, vwgt,   &
                   numflag, options,perm,   iperm)  &
      result(ierr) bind(C,name="METIS_NodeWND")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,numflag
      type(c_ptr),value           ::  xadj,adjncy
      type(c_ptr),value           ::  vwgt
      type(c_ptr),value           ::  options
      type(c_ptr),value           ::  perm
      type(c_ptr),value           ::  iperm
    end function METIS_NodeWND_f

!<  Auxillary rountines
    function METIS_EstimateMemory_f                   &
                  (nvtxs,   xadj,   adjncy, numflag,  &
                   optype,  nbytes)                   &
      result(ierr) bind(C,name="METIS_EstimateMemory")
      import
      integer(idx_t)              ::  ierr
      integer(idx_t),intent(in)   ::  nvtxs,numflag
      type(c_ptr),value           ::  xadj,adjncy
      integer(idx_t),intent(in)   ::  optype
      integer(idx_t),intent(out)  ::  nbytes
    end function METIS_EstimateMemory_f    

  end interface
  
  
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
!    if(allocated(mts%vwgt))    deallocate(mts%vwgt,stat=istat)
!    if(allocated(mts%vsize))   deallocate(mts%vsize,stat=istat)
!    if(allocated(mts%adjwgt))  deallocate(mts%adjwgt,stat=istat)
!    if(allocated(mts%tpwgts))  deallocate(mts%tpwgts,stat=istat)
!    if(allocated(mts%ubvec))   deallocate(mts%ubvec,stat=istat)
    
    if(associated(mts%options)) deallocate(mts%options,stat=istat)
    if(associated(mts%part))    deallocate(mts%part,stat=istat)
!
!    if(associated(mts%eptr))  deallocate(mts%eptr,stat=istat)
!    if(associated(mts%eind))  deallocate(mts%eind,stat=istat)
!    if(associated(mts%epart)) deallocate(mts%epart,stat=istat)
!    if(associated(mts%npart)) deallocate(mts%npart,stat=istat)
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
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  edgecut
  integer(idx_t),pointer      ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_PartGraphRecursive_f                       &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,  &
               nparts,    options_ptr,edgecut,    part_ptr)
end function METIS_PartGraphRecursive

function METIS_PartGraphKway                        &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  options,edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  edgecut
  integer(idx_t),pointer      ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_PartGraphKway_f                            &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,  &
               nparts,    options_ptr,edgecut,    part_ptr)
end function METIS_PartGraphKway

function METIS_PartGraphVKway                       &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   vsize,  wgtflag,  numflag,  &
               nparts,  options,volume,   part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  volume
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),vsize(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  volume
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,vsize_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(vsize)) then
      vsize_ptr  = c_loc(vsize)
    else
      vsize_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_PartGraphVKway_f                           &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  vsize_ptr,  wgtflag,    numflag,  &
               nparts,    options_ptr,volume,     part_ptr)
end function METIS_PartGraphVKway

function METIS_mCPartGraphRecursive                 &
              (nvtxs,   ncon,   xadj,     adjncy,   &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  options,edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_mCPartGraphRecursive_f                       &
              (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr, &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,    &
               nparts,    options_ptr,edgecut,    part_ptr)
end function METIS_mCPartGraphRecursive

function METIS_mCPartGraphKway                      &
              (nvtxs,   ncon,   xadj,     adjncy,   &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  ubvec,  options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  ubvec(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  ubvec(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  ubvec_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(ubvec)) then
      ubvec_ptr  = c_loc(ubvec)
    else
      ubvec_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_mCPartGraphKway_f                            &
              (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr, &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,    &
               nparts,    ubvec_ptr,  options_ptr,edgecut,    part_ptr)
end function METIS_mCPartGraphKway

function METIS_WPartGraphRecursive                  &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  tpwgts_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(tpwgts)) then
      tpwgts_ptr  = c_loc(tpwgts)
    else
      tpwgts_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_WPartGraphRecursive_f                      &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,  &
               nparts,    tpwgts_ptr, options_ptr,edgecut,    part_ptr)
end function METIS_WPartGraphRecursive

function METIS_WPartGraphKway                       &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  edgecut,  part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  edgecut
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  edgecut
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  tpwgts_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(tpwgts)) then
      tpwgts_ptr  = c_loc(tpwgts)
    else
      tpwgts_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_WPartGraphKway_f                           &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,  &
               nparts,    tpwgts_ptr, options_ptr,edgecut,    part_ptr)
end function METIS_WPartGraphKway

function METIS_WPartGraphVKway                      &
              (nvtxs,   xadj,   adjncy,             &
                vwgt,   adjwgt, wgtflag,  numflag,  &
               nparts,  tpwgts, options,  volume,   part) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,nparts
  integer(idx_t),intent(in)         ::  wgtflag,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  volume
!  integer(idx_t),intent(out),target ::  part(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),intent(out)        ::  volume
  integer(idx_t),pointer            ::  part(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr,adjwgt_ptr
  type(c_ptr)   ::  tpwgts_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  part_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(tpwgts)) then
      tpwgts_ptr  = c_loc(tpwgts)
    else
      tpwgts_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    
    ierr = METIS_WPartGraphVKway_f                          &
              (nvtxs,     xadj_ptr,   adjncy_ptr,           &
               vwgt_ptr,  adjwgt_ptr, wgtflag,    numflag,  &
               nparts,    tpwgts_ptr, options_ptr,volume,   part_ptr)
end function METIS_WPartGraphVKway

function METIS_EdgeND                             &
              (nvtxs,   xadj,   adjncy, numflag,  &
               options, perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out),target ::  perm(:)
!  integer(idx_t),intent(out),target ::  iperm(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  perm_ptr
  type(c_ptr)   ::  iperm_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(perm)) then
      perm_ptr  = c_loc(perm)
    else
      perm_ptr  = c_null_ptr
    endif
    if(associated(iperm)) then
      iperm_ptr  = c_loc(iperm)
    else
      iperm_ptr  = c_null_ptr
    endif
!
    ierr = METIS_EdgeND_f(nvtxs,xadj_ptr,adjncy_ptr,numflag,   &
                          options_ptr,perm_ptr,iperm_ptr)
end function METIS_EdgeND

function METIS_NodeND                             &
              (nvtxs,   xadj,   adjncy, numflag,  &
               options, perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out),target ::  perm(:)
!  integer(idx_t),intent(out),target ::  iperm(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  perm_ptr
  type(c_ptr)   ::  iperm_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(perm)) then
      perm_ptr  = c_loc(perm)
    else
      perm_ptr  = c_null_ptr
    endif
    if(associated(iperm)) then
      iperm_ptr  = c_loc(iperm)
    else
      iperm_ptr  = c_null_ptr
    endif
!
    ierr = METIS_NodeND_f(nvtxs,xadj_ptr,adjncy_ptr,numflag,   &
                          options_ptr,perm_ptr,iperm_ptr)
end function METIS_NodeND

function METIS_NodeWND                          &
              (nvtxs,   xadj,   adjncy, vwgt,   &
               numflag, options,perm,   iperm) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out),target ::  perm(:)
!  integer(idx_t),intent(out),target ::  iperm(:)
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:)
  integer(idx_t),pointer            ::  options(:)
  integer(idx_t),pointer            ::  perm(:)
  integer(idx_t),pointer            ::  iperm(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
  type(c_ptr)   ::  vwgt_ptr
  type(c_ptr)   ::  options_ptr
  type(c_ptr)   ::  perm_ptr
  type(c_ptr)   ::  iperm_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(perm)) then
      perm_ptr  = c_loc(perm)
    else
      perm_ptr  = c_null_ptr
    endif
    if(associated(iperm)) then
      iperm_ptr  = c_loc(iperm)
    else
      iperm_ptr  = c_null_ptr
    endif
!
    ierr = METIS_NodeWND_f(nvtxs,xadj_ptr,adjncy_ptr,vwgt_ptr,   &
                           numflag,options_ptr,perm_ptr,iperm_ptr)
end function METIS_NodeWND

function METIS_EstimateMemory                      &
              (nvtxs,   xadj,   adjncy, numflag,  &
               optype,  nbytes) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,numflag,optype
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
  integer(idx_t),intent(out)        ::  nbytes
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
!
  type(c_ptr)   ::  xadj_ptr,adjncy_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
!
    ierr = METIS_EstimateMemory_f(nvtxs,xadj_ptr,adjncy_ptr,numflag,   &
                          optype,nbytes)
end function METIS_EstimateMemory

!<  Extend METIS_PartGraphRecursive with super vertex information
function MetisF_SuperVertexGraph                      &
                  (nvtxs,   ncon,   xadj,   adjncy,   &
                    vwgt,   vsize,  adjwgt, tpwgts,   &
                   ubvec,   nparts,                   &
                   nsvtxs,  svptr,  svind,  mts) result(ierr)
  integer(idx_t),intent(in)         ::  nvtxs,ncon
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:),ubvec(:)
!  integer(idx_t),intent(in)         ::  nparts
!  integer(idx_t),intent(in)         ::  nsvtxs
!  integer(idx_t),intent(in),target  ::  svptr(:),svind(:)
  type(METIS4),intent(inout)        ::  mts
  integer(idx_t)                    ::  ierr
!
  integer(idx_t),pointer            ::  xadj(:),adjncy(:)
  integer(idx_t),pointer            ::  vwgt(:),vsize(:),adjwgt(:)
  real(real_t),pointer              ::  tpwgts(:),ubvec(:)
  integer(idx_t),intent(in)         ::  nparts
  integer(idx_t),intent(in)         ::  nsvtxs
  integer(idx_t),pointer            ::  svptr(:),svind(:)
!
  integer(idx_t)  ::  i,ii,j,k,istat,nvtxsLeft,count,maxNum,vs
  integer(idx_t),allocatable  ::  help(:),mark_super(:),mark_left(:)
  integer(idx_t),allocatable  ::  svptr_bc(:),svind_bc(:)
  type(CONNECT),allocatable   ::  vtxNew(:)
!
  type(c_ptr)         ::  ptr
  integer(c_intptr_t) ::  c_ptr_adr
!
    if(.not.associated(xadj).or.    &
       .not.associated(adjncy).or.  &
       .not.associated(svptr).or.   &
       .not.associated(svind)) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
!
    ierr = METIS_OK
!   Check super vertex validality
    allocate(help(nvtxs),stat=istat)
    help(:) = 0
    do i=1,nsvtxs
      do j=svptr(i),svptr(i+1)-1
        if(help(svind(j)) /= 0) then
          ierr = METIS_ERROR_SUPERVERTEX
          return      
        else
          help(svind(j)) = i
        endif
      enddo
    enddo
!   Vertex index mapping old to new & temp memory to save new graph
    mts%nvtxs = nvtxs + nsvtxs - (svptr(nsvtxs+1)-1)
    allocate(vtxNew(mts%nvtxs),stat=istat)
    nvtxsLeft = nvtxs - (svptr(nsvtxs+1)-1)
    allocate(mts%idxOld2New(nvtxs),stat=istat)
    mts%idxOld2New(:) = 0
    count = 0
    do i=1,nvtxs
      if(help(i) == 0) then
        count = count + 1
        mts%idxOld2New(i) = count
        maxNum = xadj(i+1) - xadj(i)
        allocate(vtxNew(count)%idx(maxNum),stat=istat)
      else
        mts%idxOld2New(i) = nvtxsLeft + help(i)
        maxNum = 0
        do j=xadj(i),xadj(i+1)-1
          if(help(adjncy(j)) /= help(i)) then
            vtxNew(nvtxsLeft + help(i))%num = vtxNew(nvtxsLeft + help(i))%num + 1
          endif
        enddo
      endif
    enddo
    do i=1,nsvtxs
      allocate(vtxNew(nvtxsLeft+i)%idx(vtxNew(nvtxsLeft+i)%num),stat=istat)
      vtxNew(nvtxsLeft+i)%num = 0
    enddo
!   Get boundary vertice index of super vertex, which would probably communicate with other newly generated graph vertices
    allocate(svptr_bc(nsvtxs+1),stat=istat)
    svptr_bc(:) = 0
    do i=1,nsvtxs
      next_vertex_in_super_vertex1 : do j=svptr(i),svptr(i+1)-1
        do k=xadj(svind(j)),xadj(svind(j)+1)-1
          if(help(svind(j)) /= help(adjncy(k))) then
            svptr_bc(i+1) = svptr_bc(i+1) + 1
            cycle next_vertex_in_super_vertex1
          endif
        enddo
      enddo next_vertex_in_super_vertex1
    enddo
    svptr_bc(1) = 1
    do i=1,nsvtxs
      svptr_bc(i+1) = svptr_bc(i+1) + svptr_bc(i)
    enddo
    allocate(svind_bc(svptr_bc(nsvtxs+1)-1),stat=istat)
    count = 0
    do i=1,nsvtxs
      next_vertex_in_super_vertex2 : do j=svptr(i),svptr(i+1)-1
        do k=xadj(svind(j)),xadj(svind(j)+1)-1
          if(help(svind(j)) /= help(adjncy(k))) then
            count = count + 1
            svind_bc(count) = svind(j)
            cycle next_vertex_in_super_vertex2
          endif
        enddo 
      enddo next_vertex_in_super_vertex2
    enddo     
!
    allocate(mark_super(nsvtxs),mark_left(nvtxs),stat=istat)
    count = 0
    do i=1,nvtxs
      mark_left(:) = 0
      if(help(i) == 0) then
        mark_super(:) = 0
        count = count + 1     
        do j=xadj(i),xadj(i+1)-1
          if(help(adjncy(j)) == 0) then
            vtxNew(count)%num = vtxNew(count)%num + 1
            vtxNew(count)%idx(vtxNew(count)%num) = mts%idxOld2New(adjncy(j))
          else
            if(mark_super(help(adjncy(j))) == 0) then
              mark_super(help(adjncy(j))) = i
              vtxNew(count)%num = vtxNew(count)%num + 1
              vtxNew(count)%idx(vtxNew(count)%num) = nvtxsLeft + help(adjncy(j))
            endif
          endif
        enddo
      endif
    enddo
    
    do ii=1,nsvtxs
      mark_super(:) = 0
      mark_left(:) = 0
      do i=svptr(ii),svptr(ii+1)-1
        do j=xadj(svind(i)),xadj(svind(i)+1)-1
          if(help(adjncy(j)) == help(svind(i))) cycle
          if(help(adjncy(j)) == 0) then
            if(mark_left(adjncy(j)) == 0) then
              mark_left(adjncy(j)) = 1
              vtxNew(nvtxsLeft + help(svind(i)))%num = vtxNew(nvtxsLeft + help(svind(i)))%num + 1
              vtxNew(nvtxsLeft + help(svind(i)))%idx(vtxNew(nvtxsLeft + help(svind(i)))%num) = mts%idxOld2New(adjncy(j))
            endif
          else
            if(mark_super(help(adjncy(j))) == 0) then
              mark_super(help(adjncy(j))) = 1
              vtxNew(nvtxsLeft + help(svind(i)))%num = vtxNew(nvtxsLeft + help(svind(i)))%num + 1
              vtxNew(nvtxsLeft + help(svind(i)))%idx(vtxNew(nvtxsLeft + help(svind(i)))%num) = nvtxsLeft + help(adjncy(j))
            endif
          endif
        enddo
      enddo
    enddo
!    deallocate(mark_super,stat=istat)
!    deallocate(mark_left,stat=istat)
!   New graph xadj(:)
    allocate(mts%xadj(mts%nvtxs+1),stat=istat)
    mts%xadj(1) = 1
    do i=1,mts%nvtxs
      mts%xadj(i+1) = mts%xadj(i) + vtxNew(i)%num
    enddo
!   New graph adjncy(:)
    allocate(mts%adjncy(mts%xadj(mts%nvtxs+1)-1),stat=istat)
    do i=1,mts%nvtxs
      do j=1,vtxNew(i)%num
        mts%adjncy(mts%xadj(i)+j-1) = vtxNew(i)%idx(j)
      enddo
    enddo
!   New graph vwgt(:)
    ptr = c_loc(vwgt)  
    c_ptr_adr = transfer(ptr,c_ptr_adr)
    if(.not.associated(vwgt)) then
!    if(size(vwgt) == 0) then
      mts%ncon = 1
    else
      mts%ncon = ncon + 1
    endif
    allocate(mts%vwgt(mts%ncon*mts%nvtxs),stat=istat)
    mts%vwgt(:) = 0
    do i=1,nvtxs
      if(associated(vwgt)) then
!      if(size(vwgt) > 0) then
        do j=1,ncon
          mts%vwgt((mts%idxOld2New(i)-1)*mts%ncon + j) = mts%vwgt((mts%idxOld2New(i)-1)*mts%ncon + j) + vwgt((i-1)*ncon+j)      
        enddo
      else
        mts%vwgt(mts%idxOld2New(i)) = mts%vwgt(mts%idxOld2New(i)) + 1
      endif   
    enddo
!   New graph vsize(:) must be provided since the data amount of the newly condensed super vertice have changes
    allocate(mts%vsize(mts%nvtxs),stat=istat)
    mts%vsize(:) = 0   
    do i=1,nvtxs
      if(help(i) /= 0) cycle
      vs = 1
      if(associated(vsize)) vs = vsize(i)
!      if(size(vsize) > 0) vs = vsize(i)
      mts%vsize(mts%idxOld2New(i)) = vs
    enddo
    do i=1,nsvtxs
      do j=svptr_bc(i),svptr_bc(i+1)-1
        vs = 1
        if(associated(vsize)) vs = vsize(svind_bc(j))
!        if(size(vsize) > 0) vs = vsize(svind_bc(j))
        mts%vsize(nvtxsLeft+i) = mts%vsize(nvtxsLeft+i) + vs
      enddo
    enddo
    deallocate(svptr_bc,stat=istat)
    deallocate(svind_bc,stat=istat)
!   New graph adjwgt(:) if necessary
    if(associated(adjwgt)) then
!    if(size(adjwgt) > 0) then
      allocate(mts%adjwgt(size(mts%adjncy)),stat=istat)
!     Non-condensed original vertice
      do i=1,nvtxs
        if(help(i) == 0) then
          mark_super(:) = 0
          do j=xadj(i),xadj(i+1)-1
            if(help(adjncy(j)) /= 0)  &
            mark_super(help(adjncy(j))) = mark_super(help(adjncy(j))) + adjwgt(j)
          enddo
          count = 0
          do j=xadj(i),xadj(i+1)-1
            if(help(adjncy(j)) == 0) then
              count = count + 1
              mts%adjwgt(mts%xadj(mts%idxOld2New(i)) - 1 + count) = adjwgt(j)
            else
              if(mark_super(help(adjncy(j))) > 0) then
                count = count + 1
                mts%adjwgt(mts%xadj(mts%idxOld2New(i)) - 1 + count) = mark_super(help(adjncy(j)))
                mark_super(help(adjncy(j))) = 0
              endif
            endif
          enddo
        endif
      enddo
!     Condensed super vertice
      do ii=1,nsvtxs
        mark_super(:) = 0
        mark_left(:) = 0
        do i=svptr(ii),svptr(ii+1)-1
          do j=xadj(svind(i)),xadj(svind(i)+1)-1
            if(help(adjncy(j)) == help(svind(i))) cycle
            if(help(adjncy(j)) == 0) then
              mark_left(adjncy(j)) = mark_left(adjncy(j)) + adjwgt(j)
            else
              mark_super(help(adjncy(j))) = mark_super(help(adjncy(j))) + adjwgt(j)
            endif
          enddo
        enddo
!  
        count = 0
        do i=svptr(ii),svptr(ii+1)-1
          do j=xadj(svind(i)),xadj(svind(i)+1)-1
            if(help(adjncy(j)) == help(svind(i))) cycle
            if(help(adjncy(j)) == 0) then
              if(mark_left(adjncy(j)) > 0) then
                count = count + 1
                mts%adjwgt(mts%xadj(nvtxsLeft+ii)+count-1) = mark_left(adjncy(j))
                mark_left(adjncy(j)) = 0
              endif
            elseif(help(adjncy(j)) /= 0) then
              if(mark_super(help(adjncy(j))) > 0) then
                count = count + 1
                mts%adjwgt(mts%xadj(nvtxsLeft+ii)+count-1) = mark_super(help(adjncy(j)))
                mark_super(help(adjncy(j))) = 0
              endif
            else
              continue
            endif
          enddo
        enddo      
      enddo
    endif
    deallocate(mark_super,stat=istat)
    deallocate(mark_left,stat=istat)
!   New graph tpwgts(:) if necessary
    if(mts%ncon > 1.and.associated(tpwgts)) then
      if(size(tpwgts) == nparts*ncon) then
        allocate(mts%tpwgts(nparts*mts%ncon),stat=istat)
        do i=1,nparts
          do j=1,mts%ncon
            if(j <= ncon) then
              mts%tpwgts((i-1)*mts%ncon + j) = tpwgts((i-1)*ncon + j)
            else
              mts%tpwgts((i-1)*mts%ncon + j) = 1.0D0/nparts
            endif
          enddo
        enddo
      endif
    endif   
!   New graph ubvec(:) if necessary
    ptr = c_loc(ubvec)  
    c_ptr_adr = transfer(ptr,c_ptr_adr)
    if(mts%ncon > 1.and.associated(ubvec)) then
      allocate(mts%ubvec(mts%ncon),stat=istat)
      do i=1,mts%ncon
        if(i <= size(ubvec)) then
          mts%ubvec(i) = ubvec(i)
        else
          mts%ubvec(i) = 1.01D0
        endif
      enddo   
    endif
!
    deallocate(help,stat=istat)
    call freeConnect(vtxNew)
    deallocate(vtxNew,stat=istat)
end function MetisF_SuperVertexGraph

function METIS_PartGraphRecursive_SuperVertex               &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                    vwgt,   vsize,  adjwgt, nparts, tpwgts, &
                   ubvec,   options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:),ubvec(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  objval
!  integer(idx_t),intent(out),target ::  part(:)
!  integer(idx_t),intent(in)         ::  numflag
!  integer(idx_t),intent(in)         ::  nsvtxs
!  integer(idx_t),intent(in),target  ::  svptr(:),svind(:)
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),vsize(:),adjwgt(:)
  real(real_t),pointer        ::  tpwgts(:),ubvec(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  objval
  integer(idx_t),pointer      ::  part(:)
  integer(idx_t),intent(in)   ::  numflag
  integer(idx_t),intent(in)   ::  nsvtxs
  integer(idx_t),pointer      ::  svptr(:),svind(:)
!  
  type(c_ptr)     ::  xadj_ptr,adjncy_ptr
  type(c_ptr)     ::  vwgt_ptr,vsize_ptr,adjwgt_ptr
  type(c_ptr)     ::  tpwgts_ptr,ubvec_ptr
  type(c_ptr)     ::  options_ptr
  type(c_ptr)     ::  part_ptr
  type(c_ptr)     ::  svptr_ptr,svind_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(vsize)) then
      vsize_ptr  = c_loc(vsize)
    else
      vsize_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(tpwgts)) then
      tpwgts_ptr  = c_loc(tpwgts)
    else
      tpwgts_ptr  = c_null_ptr
    endif
    if(associated(ubvec)) then
      ubvec_ptr  = c_loc(ubvec)
    else
      ubvec_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    if(associated(svptr)) then
      svptr_ptr  = c_loc(svptr)
    else
      svptr_ptr  = c_null_ptr
    endif
    if(associated(svind)) then
      svind_ptr  = c_loc(svind)
    else
      svind_ptr  = c_null_ptr
    endif
!
    ierr = METIS_PartGraphRecursive_SuperVertex_f                         &
                  (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr,         &
                   vwgt_ptr,  vsize_ptr,  adjwgt_ptr, nparts, tpwgts_ptr, &
                   ubvec_ptr, options_ptr,objval,     part_ptr,           &
                   numflag,                                               &
                   nsvtxs,    svptr_ptr,  svind_ptr)
end function METIS_PartGraphRecursive_SuperVertex

function METIS_PartGraphRecursive_SuperVertex_pgi           &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                                            nparts,         &
                            options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:),ubvec(:)
  integer(idx_t),intent(in),target  ::  options(:)
  integer(idx_t),intent(out)        ::  objval
  integer(idx_t),intent(out),target ::  part(:)
  integer(idx_t),intent(in)         ::  numflag
  integer(idx_t),intent(in)         ::  nsvtxs
  integer(idx_t),intent(in),target  ::  svptr(:),svind(:)
!
  type(c_ptr)     ::  xadj_ptr,adjncy_ptr
  type(c_ptr)     ::  vwgt_ptr,vsize_ptr,adjwgt_ptr
  type(c_ptr)     ::  tpwgts_ptr,ubvec_ptr
  type(c_ptr)     ::  options_ptr
  type(c_ptr)     ::  part_ptr
  type(c_ptr)     ::  svptr_ptr,svind_ptr
!
    xadj_ptr    = c_loc(xadj)
    adjncy_ptr  = c_loc(adjncy)
!    vwgt_ptr    = c_loc(vwgt)
!    vsize_ptr   = c_loc(vsize)
!    adjwgt_ptr  = c_loc(adjwgt)
!    tpwgts_ptr  = c_loc(tpwgts)
!    ubvec_ptr   = c_loc(ubvec)
    vwgt_ptr    = c_null_ptr
    vsize_ptr   = c_null_ptr
    adjwgt_ptr  = c_null_ptr
    tpwgts_ptr  = c_null_ptr
    ubvec_ptr   = c_null_ptr
    options_ptr = c_loc(options)
    part_ptr    = c_loc(part)
    svptr_ptr   = c_loc(svptr)
    svind_ptr   = c_loc(svind)
!
    ierr = METIS_PartGraphRecursive_SuperVertex_f                         &
                  (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr,         &
                   vwgt_ptr,  vsize_ptr,  adjwgt_ptr, nparts, tpwgts_ptr, &
                   ubvec_ptr, options_ptr,objval,     part_ptr,           &
                   numflag,                                               &
                   nsvtxs,    svptr_ptr,  svind_ptr)
end function METIS_PartGraphRecursive_SuperVertex_pgi

!<  Interface for calling from c
function METIS_PartGraphRecursive_SuperVertex_f             &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                    vwgt,   vsize,  adjwgt, nparts, tpwgts, &
                   ubvec,   options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)     &
                   bind(C,name="METIS_PartGraphRecursive_SuperVertex")
  integer(idx_t)              ::  ierr
  integer(idx_t),intent(in)   ::  nvtxs,ncon,nparts
  type(c_ptr),value           ::  xadj,adjncy
  type(c_ptr),value           ::  vwgt,vsize,adjwgt
  type(c_ptr),value           ::  tpwgts,ubvec
  type(c_ptr),value           ::  options
  integer(idx_t),intent(out)  ::  objval
  type(c_ptr),value           ::  part
  integer(idx_t),intent(in)   ::  numflag
  integer(idx_t),intent(in)   ::  nsvtxs
  type(c_ptr),value           ::  svptr,svind

!
  integer(idx_t),pointer      ::  xadj_p(:)=>NULL(),adjncy_p(:)=>NULL()
  integer(idx_t),pointer      ::  vwgt_p(:)=>NULL(),vsize_p(:)=>NULL(),adjwgt_p(:)=>NULL()
  real(real_t),pointer        ::  tpwgts_p(:)=>NULL(),ubvec_p(:)=>NULL()
  integer(idx_t),pointer      ::  options_p(:)=>NULL()
  integer(idx_t),pointer      ::  part_p(:)=>NULL()
  integer(idx_t),pointer      ::  svptr_p(:)=>NULL(),svind_p(:)=>NULL()
  integer(c_intptr_t)         ::  c_ptr_adr
  integer(idx_t)::  i,istat
  type(METIS4)  ::  tmp
!
!   Check C_PTR address to see if they are allocated
    c_ptr_adr = transfer(xadj,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(xadj,xadj_p,shape=[nvtxs+1])
    c_ptr_adr = transfer(adjncy,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(adjncy,adjncy_p,shape=[xadj_p(nvtxs+1)-1])
!
    c_ptr_adr = transfer(vwgt,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(vwgt,vwgt_p,shape=[ncon*nvtxs])
    c_ptr_adr = transfer(vsize,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(vsize,vsize_p,shape=[nvtxs])
    c_ptr_adr = transfer(adjwgt,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(adjwgt,adjwgt_p,shape=[xadj_p(nvtxs+1)-1])
!
    c_ptr_adr = transfer(tpwgts,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(tpwgts,tpwgts_p,shape=[ncon*nparts])
    c_ptr_adr = transfer(ubvec,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(ubvec,ubvec_p,shape=[ncon])
!
    c_ptr_adr = transfer(options,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(options,options_p,shape=[5])

    c_ptr_adr = transfer(part,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(part,part_p,shape=[nvtxs])

    c_ptr_adr = transfer(svptr,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(svptr,svptr_p,shape=[nsvtxs+1])

    c_ptr_adr = transfer(svind,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(svind,svind_p,shape=[svptr_p(nsvtxs+1)-1])
!    
    call MetisF_Init(tmp)
    ierr = MetisF_SuperVertexGraph                        &
                  (nvtxs,   ncon,   xadj_p,   adjncy_p,   &
                   vwgt_p,  vsize_p,adjwgt_p, tpwgts_p,   &
                   ubvec_p, nparts,                       &
                   nsvtxs,  svptr_p,  svind_p,  tmp)
!    if(ierr /= METIS_OK) return
        
    allocate(tmp%part(tmp%nvtxs),stat=istat)
!   
    call MetisF_PartGraphRecursive_Type(tmp)
    select case(tmp%optype)
    case(METIS_TYPE_PartGraphRecursive)
      ierr = METIS_PartGraphRecursive                           &
                  (tmp%nvtxs, tmp%xadj,   tmp%adjncy,           &
                   tmp%vwgt,  tmp%adjwgt, tmp%wgtflag,numflag,  &
                   nparts,    options_p,  objval,     tmp%part)    
    case(METIS_TYPE_mCPartGraphRecursive)
      ierr = METIS_mCPartGraphRecursive                           &
                  (tmp%nvtxs, tmp%ncon,   tmp%xadj,   tmp%adjncy, &
                   tmp%vwgt,  tmp%adjwgt, tmp%wgtflag,numflag,    &
                   nparts,    options_p,  objval,     tmp%part)     
    case(METIS_TYPE_WPartGraphRecursive)
      ierr = METIS_WPartGraphRecursive                            &
                  (tmp%nvtxs, tmp%xadj,   tmp%adjncy,             &
                   tmp%vwgt,  tmp%adjwgt, tmp%wgtflag,numflag,    &
                   nparts,    tmp%tpwgts, options_p,  objval,     tmp%part)    
    case default
      ierr = METIS_ERROR_SUPERVERTEX
      return    
    end select

!    if(ierr /= METIS_OK) return
    if(nparts == 1) then
      part_p(:) = 1
    else
!$omp parallel shared(nvtxs,part_p,tmp) private(i)
!$omp do
    do i=1,nvtxs
      part_p(i) = tmp%part(tmp%idxOld2New(i))
    enddo
!$omp end do
!$omp end parallel
    endif
    
    call MetisF_Free(tmp)
end function METIS_PartGraphRecursive_SuperVertex_f

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

function METIS_PartGraphKway_SuperVertex                    &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                    vwgt,   vsize,  adjwgt, nparts, tpwgts, &
                   ubvec,   options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
!  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:),ubvec(:)
!  integer(idx_t),intent(in),target  ::  options(:)
!  integer(idx_t),intent(out)        ::  objval
!  integer(idx_t),intent(out),target ::  part(:)
!  integer(idx_t),intent(in)         ::  numflag
!  integer(idx_t),intent(in)         ::  nsvtxs
!  integer(idx_t),intent(in),target  ::  svptr(:),svind(:)
!
  integer(idx_t),pointer      ::  xadj(:),adjncy(:)
  integer(idx_t),pointer      ::  vwgt(:),vsize(:),adjwgt(:)
  real(real_t),pointer        ::  tpwgts(:),ubvec(:)
  integer(idx_t),pointer      ::  options(:)
  integer(idx_t),intent(out)  ::  objval
  integer(idx_t),pointer      ::  part(:)
  integer(idx_t),intent(in)   ::  numflag
  integer(idx_t),intent(in)   ::  nsvtxs
  integer(idx_t),pointer      ::  svptr(:),svind(:)
!  
  type(c_ptr)     ::  xadj_ptr,adjncy_ptr
  type(c_ptr)     ::  vwgt_ptr,vsize_ptr,adjwgt_ptr
  type(c_ptr)     ::  tpwgts_ptr,ubvec_ptr
  type(c_ptr)     ::  options_ptr
  type(c_ptr)     ::  part_ptr
  type(c_ptr)     ::  svptr_ptr,svind_ptr
!
    if(associated(xadj)) then
      xadj_ptr  = c_loc(xadj)
    else
      xadj_ptr  = c_null_ptr
    endif
    if(associated(adjncy)) then
      adjncy_ptr  = c_loc(adjncy)
    else
      adjncy_ptr  = c_null_ptr
    endif
    if(associated(vwgt)) then
      vwgt_ptr  = c_loc(vwgt)
    else
      vwgt_ptr  = c_null_ptr
    endif
    if(associated(vsize)) then
      vsize_ptr  = c_loc(vsize)
    else
      vsize_ptr  = c_null_ptr
    endif
    if(associated(adjwgt)) then
      adjwgt_ptr  = c_loc(adjwgt)
    else
      adjwgt_ptr  = c_null_ptr
    endif
    if(associated(tpwgts)) then
      tpwgts_ptr  = c_loc(tpwgts)
    else
      tpwgts_ptr  = c_null_ptr
    endif
    if(associated(ubvec)) then
      ubvec_ptr  = c_loc(ubvec)
    else
      ubvec_ptr  = c_null_ptr
    endif
    if(associated(options)) then
      options_ptr  = c_loc(options)
    else
      options_ptr  = c_null_ptr
    endif
    if(associated(part)) then
      part_ptr  = c_loc(part)
    else
      part_ptr  = c_null_ptr
    endif
    if(associated(svptr)) then
      svptr_ptr  = c_loc(svptr)
    else
      svptr_ptr  = c_null_ptr
    endif
    if(associated(svind)) then
      svind_ptr  = c_loc(svind)
    else
      svind_ptr  = c_null_ptr
    endif
!
    ierr = METIS_PartGraphKway_SuperVertex_f                              &
                  (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr,         &
                   vwgt_ptr,  vsize_ptr,  adjwgt_ptr, nparts, tpwgts_ptr, &
                   ubvec_ptr, options_ptr,objval,     part_ptr,           &
                   numflag,                                               &
                   nsvtxs,    svptr_ptr,  svind_ptr)
end function METIS_PartGraphKway_SuperVertex

function METIS_PartGraphKway_SuperVertex_pgi                &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                                            nparts,         &
                            options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)
  integer(idx_t)                    ::  ierr
  integer(idx_t),intent(in)         ::  nvtxs,ncon,nparts
  integer(idx_t),intent(in),target  ::  xadj(:),adjncy(:)
!  integer(idx_t),intent(in),target  ::  vwgt(:),vsize(:),adjwgt(:)
!  real(real_t),intent(in),target    ::  tpwgts(:),ubvec(:)
  integer(idx_t),intent(in),target  ::  options(:)
  integer(idx_t),intent(out)        ::  objval
  integer(idx_t),intent(out),target ::  part(:)
  integer(idx_t),intent(in)         ::  numflag
  integer(idx_t),intent(in)         ::  nsvtxs
  integer(idx_t),intent(in),target  ::  svptr(:),svind(:)
!
  type(c_ptr)     ::  xadj_ptr,adjncy_ptr
  type(c_ptr)     ::  vwgt_ptr,vsize_ptr,adjwgt_ptr
  type(c_ptr)     ::  tpwgts_ptr,ubvec_ptr
  type(c_ptr)     ::  options_ptr
  type(c_ptr)     ::  part_ptr
  type(c_ptr)     ::  svptr_ptr,svind_ptr
!
    xadj_ptr    = c_loc(xadj)
    adjncy_ptr  = c_loc(adjncy)
!    vwgt_ptr    = c_loc(vwgt)
!    vsize_ptr   = c_loc(vsize)
!    adjwgt_ptr  = c_loc(adjwgt)
!    tpwgts_ptr  = c_loc(tpwgts)
!    ubvec_ptr   = c_loc(ubvec)
    vwgt_ptr    = c_null_ptr
    vsize_ptr   = c_null_ptr
    adjwgt_ptr  = c_null_ptr
    tpwgts_ptr  = c_null_ptr
    ubvec_ptr   = c_null_ptr
    options_ptr = c_loc(options)
    part_ptr    = c_loc(part)
    svptr_ptr   = c_loc(svptr)
    svind_ptr   = c_loc(svind)
!
    ierr = METIS_PartGraphKway_SuperVertex_f                              &
                  (nvtxs,     ncon,       xadj_ptr,   adjncy_ptr,         &
                   vwgt_ptr,  vsize_ptr,  adjwgt_ptr, nparts, tpwgts_ptr, &
                   ubvec_ptr, options_ptr,objval,     part_ptr,           &
                   numflag,                                               &
                   nsvtxs,    svptr_ptr,  svind_ptr)
end function METIS_PartGraphKway_SuperVertex_pgi

function METIS_PartGraphKway_SuperVertex_f                  &
                  (nvtxs,   ncon,   xadj,   adjncy,         &
                    vwgt,   vsize,  adjwgt, nparts, tpwgts, &
                   ubvec,   options,objval, part,           &
                   numflag,                                 &
                   nsvtxs,  svptr,  svind) result(ierr)     &
                   bind(C,name="METIS_PartGraphKway_SuperVertex")
  integer(idx_t)              ::  ierr
  integer(idx_t),intent(in)   ::  nvtxs,ncon,nparts
  type(c_ptr),value           ::  xadj,adjncy
  type(c_ptr),value           ::  vwgt,vsize,adjwgt
  type(c_ptr),value           ::  tpwgts,ubvec
  type(c_ptr),value           ::  options
  integer(idx_t),intent(out)  ::  objval
  type(c_ptr),value           ::  part
  integer(idx_t),intent(in)   ::  numflag
  integer(idx_t),intent(in)   ::  nsvtxs
  type(c_ptr),value           ::  svptr,svind
!
 
  integer(idx_t),pointer      ::  xadj_p(:)=>NULL(),adjncy_p(:)=>NULL()
  integer(idx_t),pointer      ::  vwgt_p(:)=>NULL(),vsize_p(:)=>NULL(),adjwgt_p(:)=>NULL()
  real(real_t),pointer        ::  tpwgts_p(:)=>NULL(),ubvec_p(:)=>NULL()
  integer(idx_t),pointer      ::  options_p(:)=>NULL()
  integer(idx_t),pointer      ::  part_p(:)=>NULL()
  integer(idx_t),pointer      ::  svptr_p(:)=>NULL(),svind_p(:)=>NULL()
  integer(c_intptr_t)         ::  c_ptr_adr
  integer(idx_t)::  i,istat
  type(METIS4)  ::  tmp
!
!   Check C_PTR address to see if they are allocated
    c_ptr_adr = transfer(xadj,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(xadj,xadj_p,shape=[nvtxs+1])
    c_ptr_adr = transfer(adjncy,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(adjncy,adjncy_p,shape=[xadj_p(nvtxs+1)-1])
!
    c_ptr_adr = transfer(vwgt,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(vwgt,vwgt_p,shape=[ncon*nvtxs])
    c_ptr_adr = transfer(vsize,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(vsize,vsize_p,shape=[nvtxs])
    c_ptr_adr = transfer(adjwgt,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(adjwgt,adjwgt_p,shape=[xadj_p(nvtxs+1)-1])
!
    c_ptr_adr = transfer(tpwgts,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(tpwgts,tpwgts_p,shape=[ncon*nparts])
    c_ptr_adr = transfer(ubvec,c_ptr_adr)
    if(c_ptr_adr /= 0) call c_f_pointer(ubvec,ubvec_p,shape=[ncon])
!
    c_ptr_adr = transfer(options,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(options,options_p,shape=[5])

    c_ptr_adr = transfer(part,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(part,part_p,shape=[nvtxs])

    c_ptr_adr = transfer(svptr,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(svptr,svptr_p,shape=[nsvtxs+1])

    c_ptr_adr = transfer(svind,c_ptr_adr)
    if(c_ptr_adr == 0) then
      ierr = METIS_ERROR_SUPERVERTEX
      return
    endif
    call c_f_pointer(svind,svind_p,shape=[svptr_p(nsvtxs+1)-1])
!    
    call MetisF_Init(tmp)
    ierr = MetisF_SuperVertexGraph                        &
                  (nvtxs,   ncon,   xadj_p,   adjncy_p,   &
                   vwgt_p,  vsize_p,adjwgt_p, tpwgts_p,   &
                   ubvec_p, nparts,                       &
                   nsvtxs,  svptr_p,  svind_p,  tmp)
!    if(ierr /= METIS_OK) return
        
    allocate(tmp%part(tmp%nvtxs),stat=istat)
    call MetisF_PartGraphKway_Type(tmp)
    select case(tmp%optype)
    case(METIS_TYPE_PartGraphKway)
      ierr = METIS_PartGraphKway                                    &
                  (tmp%nvtxs,   tmp%xadj,   tmp%adjncy,             &
                   tmp%vwgt,    tmp%adjwgt, tmp%wgtflag,numflag,    &
                   nparts,      options_p,  objval,     tmp%part)
    
    case(METIS_TYPE_PartGraphVKway)
      ierr = METIS_PartGraphVKway                                   &
                  (tmp%nvtxs,   tmp%xadj,   tmp%adjncy,             &
                   tmp%vwgt,    tmp%vsize,  tmp%wgtflag,numflag,    &
                   nparts,      options_p,  objval,     tmp%part)
    
    case(METIS_TYPE_mCPartGraphKway)
      ierr = METIS_mCPartGraphKway                                  &
                  (tmp%nvtxs,   tmp%ncon,   tmp%xadj,   tmp%adjncy, &
                   tmp%vwgt,    tmp%adjwgt, tmp%wgtflag,numflag,    &
                   nparts,      tmp%ubvec,  options_p,  objval,   tmp%part)
    
    case(METIS_TYPE_WPartGraphKway)
      ierr = METIS_WPartGraphKway                                   &
                  (tmp%nvtxs,   tmp%xadj,   tmp%adjncy,             &
                   tmp%vwgt,    tmp%adjwgt, tmp%wgtflag,numflag,    &
                   nparts,      tmp%tpwgts, options_p,  objval,   tmp%part)
    
    case(METIS_TYPE_WPartGraphVKway)
      ierr = METIS_WPartGraphVKway                                  &
                  (tmp%nvtxs,   tmp%xadj,   tmp%adjncy,             &
                   tmp%vwgt,    tmp%vsize,  tmp%wgtflag,numflag,    &
                   nparts,      tmp%tpwgts, options_p,  objval,   tmp%part)

    case default
      ierr = METIS_ERROR_SUPERVERTEX
      return    
    end select
    
!    if(ierr /= METIS_OK) return

    if(nparts == 1) then
      part_p(:) = 1
    else
!$omp parallel shared(nvtxs,part_p,tmp) private(i)
!$omp do
    do i=1,nvtxs
      part_p(i) = tmp%part(tmp%idxOld2New(i))
    enddo
!$omp end do
!$omp end parallel
    endif
    
    call MetisF_Free(tmp)
end function METIS_PartGraphKway_SuperVertex_f

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