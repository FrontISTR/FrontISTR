#ifdef USE_VHCALL
subroutine stiffMatrixCallerWrapper(hecMAT, ctime, tincr, cstep, maxn_gauss)
  use vhcall_fortran
  use hecmw_util, only: kint, kreal, hecmwST_local_mesh, hecmwST_matrix, &
    hecmw_comm_get_comm, hecmw_comm_get_rank, hecmw_comm_get_size
  use m_fstr, only: fstr_solid, fstr_param
  use m_fstr_main, only: hecMESH, fstrSOLID
  implicit none
  integer(8) :: handle
  integer(8) :: func
  integer(8) :: args
  integer(8) :: ret
  type (hecmwST_matrix), intent(inout)      :: hecMAT       !< linear equation, its right side modified here
  real(kind=kreal),intent(in) :: ctime        !< current time
  real(kind=kreal),intent(in) :: tincr       !< time increment
  integer, intent(in)                  :: cstep       !< current step
  integer(kind=kint) :: n_mat
  integer(kind=kint) :: hecmw_PETOT,hecmw_rank,hecmw_comm,hecmw_group
  integer(kind=kint) :: maxn_gauss

  hecmw_PETOT=hecmw_comm_get_size()
  hecmw_rank=hecmw_comm_get_rank()
  hecmw_comm=hecmw_comm_get_comm()

  handle = fvhcall_install('libfistr.so')
  func = fvhcall_find(handle, 'stiffMatrixCalleeWrapper')
  args = fvhcall_args_alloc()

  ret = fvhcall_args_set(args, fvhcall_intent_in, 1, ctime)
  if(ret /= 0) then
    write(0,*) 'fvhcall_args_set for ctime failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 2, tincr)
  if(ret /= 0) then
    write(0,*) 'fvhcall_args_set for tincr failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 3, cstep)
  if(ret /= 0) then
    write(0,*) 'fvhcall_args_set for cstep failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 4, hecmw_comm)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for hecmw_comm failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 5, hecmw_PETOT)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for hecmw_PETOT failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 6, hecmw_rank)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for hecmw_rank failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 7, hecmw_group)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for hecmw_group failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 8, maxn_gauss)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for maxn_gauss failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_out, 9, hecMAT%D)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for D failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_out, 10, hecMAT%AL)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed AL', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_out, 11, hecMAT%AU)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for AU failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 12, size(hecMAT%D))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for size(hecMAT%D) failed', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 13, size(hecMAT%AL))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed size(hecMAT%AL)', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 14, hecMESH%node)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed hecMESH%node', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 15, size(hecMESH%node))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed size(hecMESH%node)', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 16, fstrSOLID%dunode)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed fstrSOLID%dunode', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 17, size(fstrSOLID%dunode))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed size(fstrSOLID%dunode)', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 18, fstrSOLID%unode)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed fstrSOLID%unode', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 19, size(fstrSOLID%unode))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed size(fstrSOLID%unode)', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 20, fstrSOLID%factor)
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed fstrSOLID%factor', ret
  endif
  ret = fvhcall_args_set(args, fvhcall_intent_in, 21, size(fstrSOLID%factor))
  if(ret /= 0)then
    write(0,*) 'fvhcall_args_set for failed size(fstrSOLID%factor)', ret
  endif

  ret = fvhcall_invoke_with_args(func, args)
  if(ret /= 0)then
    write(0,*) 'VH call failed',ret
  endif

end subroutine
#endif

#ifdef VHCALL_LIB
subroutine stiffMatrixCalleeWrapper(ctime, tincr, cstep, hecmw_comm, hecmw_PETOT, hecmw_rank, hecmw_group, maxn_gauss, &
    D, AL, AU, sizeD, sizeALAU, node, sizeNode, dunode, sizeDunode, unode, sizeUnode, factor, sizeFactor )
  use hecmw_util, only: kint, kreal, HECMW_FILENAME_LEN
  use m_fstr, only: fstr_param
  use m_fstr_main, only: hecMESH, hecMAT, fstrSOLID, fstrPR, fstrCPL, &
      hecmw_get_mesh, hecmw2fstr_mesh_conv, fstr_init, fstr_rcap_initialize
  use m_fstr_StiffMatrix, only: fstr_StiffMatrix
  use m_fstr_spring, only: fstr_AddSPRING
  implicit none
  real(kind=kreal),intent(in) :: ctime
  real(kind=kreal),intent(in) :: tincr
  integer, intent(in)         :: cstep
  character(len=HECMW_FILENAME_LEN)          :: name_ID
  integer(kind=kint) :: hecmw_PETOT,hecmw_rank,hecmw_comm,hecmw_group
  character(len=HECMW_FILENAME_LEN):: ctrlfile = "hecmw_ctrl.dat"
  integer(kind=kint) :: ierr
  integer(kind=kint) :: maxn_gauss
  integer(kind=kint) :: sizeD, sizeALAU, sizeNode, sizeDunode, sizeUnode, sizeFactor
  real(kind=kreal),intent(out) :: D(sizeD), AL(sizeALAU), AU(sizeALAU)
  real(kind=kreal),intent(in) :: node(sizeNode), dunode(sizeDunode), unode(sizeUnode), factor(sizeFactor)
  integer(kind=kint) :: i
  type(fstr_param)   :: fstrPARAM
  logical,save       :: is_first=.true.


  if( is_first) then
    call hecmw_comm_init_if(hecmw_comm, hecmw_PETOT, hecmw_rank, hecmw_group)
    call hecmw_ctrl_init_ex_if(ctrlfile, ierr)

    name_ID = 'fstrMSH'
    call hecmw_get_mesh( name_ID , hecMESH )
    call hecmw2fstr_mesh_conv( hecMESH )
    call fstr_init

    fstrSOLID%maxn_gauss=maxn_gauss
    call fstr_rcap_initialize( hecMESH, fstrPR, fstrCPL )
    is_first=.false.
  else
    hecMESH%node=node
    fstrSOLID%dunode=dunode
    fstrSOLID%unode=unode
    fstrSOLID%factor=factor
  endif

  call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, ctime, tincr )
  call fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

  do i=1, sizeD
    D(i)=hecMAT%D(i)
  enddo
  do i=1, sizeALAU
    AL(i)=hecMAT%AL(i)
    AU(i)=hecMAT%AU(i)
  enddo
end subroutine
#endif
