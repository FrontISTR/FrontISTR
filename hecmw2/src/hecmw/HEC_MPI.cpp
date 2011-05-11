//
//  HEC_MW.cpp
//
//  mpichのラッパー
//
//                  2010.03.12
//                  2009.08.24
//                  k.Takeda
#include "HEC_MPI.h"
#include "Logger.h"
using namespace pmw;

// construct & destruct
//
CHecMPI::CHecMPI()
{
    ;
}
CHecMPI::~CHecMPI()
{
    ;
}

// 初期化処理
//
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

    //--------------------   Message   -----------------------
    //  "MPI Initialize : PE "  ## "Number of Process"  ##
    //--------------------------------------------------------
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

//
//
int CHecMPI::Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
#else
    return 0;
#endif
//#ifdef HAVE_MPI
//    int rtc;
//
//	if( datatype == HECMW_INT ) {
//        rtc = MPI_Isend( buf, count, MPI_INT, dest, tag, comm, request );
//        if( rtc != MPI_SUCCESS ) {
//            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
//            goto error;
//        }
//	} else if( datatype == HECMW_DOUBLE ) {
//        rtc = MPI_Isend( buf, count, MPI_DOUBLE, dest, tag, comm, request );
//        if( rtc != MPI_SUCCESS ) {
//            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
//            goto error;
//        }
//	} else if( datatype == HECMW_CHAR ) {
//        rtc = MPI_Isend( buf, count, MPI_CHAR, dest, tag, comm, request );
//        if( rtc != MPI_SUCCESS ) {
//            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
//            goto error;
//        }
//	} else {
//        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
//        goto error;
//    }
//    return 0;
//error:
//    return -1;
//#else
//	return 0;
//#endif
}

//
//
int CHecMPI::Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#else
    return 0;
#endif
}

//
//
int CHecMPI::Wait(MPI_Request *request, MPI_Status *status)
{
#ifdef HAVE_MPI
    return MPI_Wait(request, status);
#else
    return 0;
#endif
}

//
//
int CHecMPI::Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
#else
    recvbuf = sendbuf;
    return 0;
#endif
}

//
//
int CHecMPI::Barrier(MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Barrier(comm);
#else
    return 0;
#endif
}

//
//
int CHecMPI::Abort(MPI_Comm comm, int errorcode)
{
#ifdef HAVE_MPI
    return MPI_Abort(comm, errorcode);
#else
    return 0;
#endif
}

//
//
int CHecMPI::Comm_rank(MPI_Comm comm, int* rank)
{
#ifdef HAVE_MPI
    return MPI_Comm_rank(comm, rank);
#else
    *rank=0;
    return 0;
#endif
}

// buf:送信データの先頭アドレス
// count:送信データの個数
// datatype:送信データの型
// dest:受信先ランク
// tag:メッセージ・タグ(送受信が同一)
//
int CHecMPI::Send(void* buf, int count, MPI_Datatype datatype, int dest,
                  int tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Send(buf, count, datatype, dest, tag, comm);
#else
    return 0;
#endif
}
// buf:受信データの先頭アドレス
// count:受信データの個数
// datatype:受信データの型
// source:送信元ランク
// tag:メッセージ・タグ(送受信が同一)
// status:受信状態
//
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
//
//
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

//
//
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





