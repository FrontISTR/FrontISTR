/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   HEC_MPI.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
using namespace pmw;
CHecMPI::CHecMPI()
{
    ;
}
CHecMPI::~CHecMPI()
{
    ;
}
void CHecMPI::Initialize(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#else
    numOfProcess = 1;
    myRank = 0;
#endif
    std::cout << " myRank => " << myRank << ",number of process => " << numOfProcess << std::endl;
}
void CHecMPI::Finalize()
{
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}
int CHecMPI::Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
#else
    return 0;
#endif
}
int CHecMPI::Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#else
    return 0;
#endif
}
int CHecMPI::Wait(MPI_Request *request, MPI_Status *status)
{
#ifdef HAVE_MPI
    return MPI_Wait(request, status);
#else
    return 0;
#endif
}
int CHecMPI::Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
#else
    return 0;
#endif
}
int CHecMPI::Barrier(MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Barrier(comm);
#else
    return 0;
#endif
}
int CHecMPI::Abort(MPI_Comm comm, int errorcode)
{
#ifdef HAVE_MPI
    return MPI_Abort(comm, errorcode);
#else
    return 0;
#endif
}
int CHecMPI::Comm_rank(MPI_Comm comm, int* rank)
{
#ifdef HAVE_MPI
    return MPI_Comm_rank(comm, rank);
#else
    return 0;
#endif
}
int CHecMPI::Send(void* buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Send(buf, count, datatype, dest, tag, comm);
#else
    return 0;
#endif
}
int CHecMPI::Recv(void* buf, int count, MPI_Datatype datatype, int source,
                      int tag, MPI_Comm comm, MPI_Status *status)
{
#ifdef HAVE_MPI
    return MPI_Recv(buf, count, datatype,  source, tag, comm, status);
#else
    return 0;
#endif
}
int CHecMPI::Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype,
                    void* recvbuf, int recvcount, MPI_Datatype recvtype,
                    int root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm);
#else
    return 0;
#endif
}
int CHecMPI::Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype,
                     void* recvbuf, int recvcnt, MPI_Datatype recvtype,
                     int root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
#else
    return 0;
#endif
}
