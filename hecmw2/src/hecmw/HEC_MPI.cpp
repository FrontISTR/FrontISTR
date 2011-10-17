/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/HEC_MPI.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "HEC_MPI.h"
#include "Logger.h"
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
    MPI_Comm_size(MPI_COMM_WORLD, &nNumOfProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    vstring vMessage; vMessage.resize(2);
    vuint   vNumber;  vNumber.resize(2);
    vMessage[0]= "MPI Initialize : PE "; vNumber[0]= myRank;
    vMessage[1]= " Number of Process ";  vNumber[1]= nNumOfProcess;
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Info, vMessage, vNumber);
#else
    nNumOfProcess = 1;
    myRank = 0;
#endif
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
    recvbuf = sendbuf;
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
    *rank=0;
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
int CHecMPI::Allgather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
#else
    recvbuf = sendbuf;
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
    recvbuf = sendbuf;
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
    recvbuf = sendbuf;
    return 0;
#endif
}
int CHecMPI::Bcast(void* buf, int cnt, MPI_Datatype type, int root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Bcast(buf, cnt, type, root, comm);
#else
    buf = buf;
    return 0;
#endif
}
