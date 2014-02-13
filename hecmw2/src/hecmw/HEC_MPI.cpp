/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/HEC_MPI.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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

    //--
    // MPI_UIINT, MPI_IINT 登録
    //--
#ifdef INT64
#ifdef MSVC
    MPI_Type_contiguous(2, MPI_UNSIGNED_LONG, &mnMPI_UIINT);//-- Win 64bit unsigned int
    MPI_Type_commit( &mnMPI_UIINT );
#else
    MPI_Type_contiguous(1, MPI_UNSIGNED_LONG, &mnMPI_UIINT);//-- Linux 64bit unsigned int
    MPI_Type_commit( &mnMPI_UIINT );
#endif
    MPI_Type_contiguous(1, MPI_LONG_LONG_INT, &mnMPI_IINT); //-- 64bit signed int
    MPI_Type_commit( &mnMPI_IINT );
#else
    MPI_Type_contiguous(1, MPI_UNSIGNED, &mnMPI_UIINT);//--32bit unsigned int
    MPI_Type_commit( &mnMPI_UIINT );

    MPI_Type_contiguous(1, MPI_INT, &mnMPI_IINT);//--------32bit signed int
    MPI_Type_commit( &mnMPI_IINT );
#endif  //INT64 End

    //--
    // MPI Initialize メッセージ出力
    //--
    string sMessage;
    stringstream ss;

    sMessage = "MPI Initialize : PE ";
    ss << myRank;
    sMessage += ss.str();
    sMessage += " Number of Process ";

    ss.str("");
    ss.clear();
    ss << nNumOfProcess;
    sMessage += ss.str();

    cout << "Info      " << setfill('_') << setw(66) << sMessage << endl;

#else
    nNumOfProcess= 1;
    myRank= 0;

    mnMPI_IINT= 0;
    mnMPI_UIINT= 0;
#endif //HAVE_MPI End

}
void CHecMPI::Finalize()
{
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
}

void CHecMPI::Comm_rank(MPI_Comm comm, int* rank)
{
#ifdef HAVE_MPI
    MPI_Comm_rank(comm, rank);
#else
    *rank=0;
#endif
}

iint CHecMPI::Isend(void* buf, iint count, MPI_Datatype datatype, iint dest, iint tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
#else
    return 0;
#endif
}
iint CHecMPI::Irecv(void* buf, iint count, MPI_Datatype datatype, iint source, iint tag, MPI_Comm comm, MPI_Request *request)
{
#ifdef HAVE_MPI
    return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#else
    return 0;
#endif
}
iint CHecMPI::Wait(MPI_Request *request, MPI_Status *status)
{
#ifdef HAVE_MPI
    return MPI_Wait(request, status);
#else
    return 0;
#endif
}
iint CHecMPI::Waitall(iint nNumOfReq, MPI_Request* request, MPI_Status* status)
{
#ifdef HAVE_MPI
    return MPI_Waitall(nNumOfReq, request, status);
#else
    return 0;
#endif
}
//--
// MPI_Op による操作を派生データ型(MPI_UIINT,MPI_IINT)にも適用
//--
iint CHecMPI::Allreduce(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
#ifdef HAVE_MPI

    iint nVal;
    // 基本型 : MPI_Allreduce
    if(datatype != mnMPI_UIINT && datatype != mnMPI_IINT) {
        nVal= MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    } else {
        ////switch(op){
        ////    case(MPI_SUM):
        ////        nVal= Allreduce_sum(sendbuf, recvbuf, count, datatype, comm);
        ////        break;
        ////    case(MPI_MAX):
        ////        nVal= Allreduce_max(sendbuf, recvbuf, count, datatype, comm);
        ////        break;
        ////    case(MPI_MIN):
        ////        nVal= Allreduce_min(sendbuf, recvbuf, count, datatype, comm);
        ////        break;
        ////    case(MPI_PROD):
        ////        nVal= Allreduce_prod(sendbuf, recvbuf, count, datatype, comm);
        ////        break;
        ////    default:
        ////    {
        ////        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        ////        pLogger->Info(Utility::LoggerMode::Error, "HecMPI::Allreduce, MPI_Op type mismatch");
        ////    }
        ////        break;
        ////}//end switch
        //
        // 特定コンパイラー対応
        //
        if(op==MPI_SUM)  nVal= Allreduce_sum(sendbuf, recvbuf, count, datatype, comm);
        if(op==MPI_MAX)  nVal= Allreduce_max(sendbuf, recvbuf, count, datatype, comm);
        if(op==MPI_MIN)  nVal= Allreduce_min(sendbuf, recvbuf, count, datatype, comm);
        if(op==MPI_PROD) nVal= Allreduce_prod(sendbuf, recvbuf, count, datatype, comm);
        if(op!=MPI_SUM && op!=MPI_MAX && op!=MPI_MIN && op!=MPI_PROD) {
            Utility::CLogger *pLogger= Utility::CLogger::Instance();
            pLogger->Info(Utility::LoggerMode::Error, "HecMPI::Allreduce, MPI_Op type mismatch");
            nVal= -1.0;
        }
    }//end if

    return nVal;
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
//--
// 派生データ型(UIINT,IINT) Allreduce : MPI_SUM
//--
iint CHecMPI::Allreduce_sum(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    iint nVal;
    //--
    // MPI_UIINT
    //--
    if(datatype==mnMPI_UIINT) {
        uiint *rbuf= (uiint*)malloc(count*nNumOfProcess*sizeof(uiint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        uiint *sumbuf= (uiint*)calloc(count, sizeof(uiint));
        // 総和
        for(int i=0; i < nNumOfProcess; i++) {
            for(iint ii=0; ii < count; ii++) {
                sumbuf[ii] += rbuf[i*count + ii];
            };
        };
        // recvbufへ代入
        uiint* buf= static_cast<uiint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = sumbuf[ii];

        free( sumbuf );
        free( rbuf );
    }
    //--
    // MPI_IINT
    //--
    if(datatype==mnMPI_IINT) {
        iint *rbuf= (iint*)malloc(count*nNumOfProcess*sizeof(iint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        iint *sumbuf= (iint*)calloc(count, sizeof(iint));
        // 総和
        for(int i=0; i < nNumOfProcess; i++) {
            for(iint ii=0; ii < count; ii++) {
                sumbuf[ii] += rbuf[i*count + ii];
            };
        };
        // recvbufへ代入
        iint* buf= static_cast<iint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = sumbuf[ii];

        free( sumbuf );
        free( rbuf );
    }

    return nVal;
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
//--
// 派生データ型(UIINT,IINT) Allreduce : MPI_MAX
//--
iint CHecMPI::Allreduce_max(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    iint nVal;
    //--
    // MPI_UIINT
    //--
    if(datatype==mnMPI_UIINT) {
        uiint *rbuf= (uiint*)malloc(count*nNumOfProcess*sizeof(uiint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        uiint *maxbuf= (uiint*)calloc(count, sizeof(uiint));
        // 最大値
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    maxbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    if(maxbuf[ii] < rbuf[i*count + ii]) maxbuf[ii]= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        uiint* buf= static_cast<uiint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = maxbuf[ii];

        free( maxbuf );
        free( rbuf );
    }
    //--
    // MPI_IINT
    //--
    if(datatype==mnMPI_IINT) {
        iint *rbuf= (iint*)malloc(count*nNumOfProcess*sizeof(iint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        iint *maxbuf= (iint*)calloc(count, sizeof(iint));
        // 最大値
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    maxbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    if(maxbuf[ii] < rbuf[i*count + ii]) maxbuf[ii]= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        iint* buf= static_cast<iint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = maxbuf[ii];

        free( maxbuf );
        free( rbuf );
    }

    return nVal;
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
//--
// 派生データ型(UIINT,IINT) Allreduce : MPI_MIN
//--
iint CHecMPI::Allreduce_min(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    iint nVal;
    //--
    // MPI_UIINT
    //--
    if(datatype==mnMPI_UIINT) {
        uiint *rbuf= (uiint*)malloc(count*nNumOfProcess*sizeof(uiint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        uiint *minbuf= (uiint*)calloc(count, sizeof(uiint));
        // 最小値
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    minbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    if(minbuf[ii] > rbuf[i*count + ii]) minbuf[ii]= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        uiint* buf= static_cast<uiint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = minbuf[ii];

        free( minbuf );
        free( rbuf );
    }
    //--
    // MPI_IINT
    //--
    if(datatype==mnMPI_IINT) {
        iint *rbuf= (iint*)malloc(count*nNumOfProcess*sizeof(iint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        iint *minbuf= (iint*)calloc(count, sizeof(iint));
        // 最小値
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    minbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    if(minbuf[ii] > rbuf[i*count + ii]) minbuf[ii]= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        iint* buf= static_cast<iint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = minbuf[ii];

        free( minbuf );
        free( rbuf );
    }

    return nVal;
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
//--
// 派生データ型(UIINT,IINT) Allreduce : MPI_PROD
//--
iint CHecMPI::Allreduce_prod(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    iint nVal;
    //--
    // MPI_UIINT
    //--
    if(datatype==mnMPI_UIINT) {
        uiint *rbuf= (uiint*)malloc(count*nNumOfProcess*sizeof(uiint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        uiint *prodbuf= (uiint*)calloc(count, sizeof(uiint));
        // 積
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    prodbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    prodbuf[ii]*= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        uiint* buf= static_cast<uiint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = prodbuf[ii];

        free( prodbuf );
        free( rbuf );
    }
    //--
    // MPI_IINT
    //--
    if(datatype==mnMPI_IINT) {
        iint *rbuf= (iint*)malloc(count*nNumOfProcess*sizeof(iint));
        nVal= MPI_Allgather(sendbuf, count, datatype, rbuf, count, datatype, comm);

        iint *prodbuf= (iint*)calloc(count, sizeof(iint));
        // 積
        for(int i=0; i < nNumOfProcess; i++) {
            if(i==0) {
                for(iint ii=0; ii < count; ii++) {
                    prodbuf[ii]= rbuf[i*count + ii];
                };
            } else {
                for(iint ii=0; ii < count; ii++) {
                    prodbuf[ii]*= rbuf[i*count +ii];
                };
            }// end if
        };
        // recvbufへ代入
        iint* buf= static_cast<iint*>(recvbuf);
        for(iint ii=0; ii < count; ii++) buf[ii] = prodbuf[ii];

        free( prodbuf );
        free( rbuf );
    }

    return nVal;
#else
    recvbuf = sendbuf;
    return 0;
#endif
}


iint CHecMPI::Barrier(MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Barrier(comm);
#else
    return 0;
#endif
}
iint CHecMPI::Abort(MPI_Comm comm, iint errorcode)
{
#ifdef HAVE_MPI
    return MPI_Abort(comm, errorcode);
#else
    return 0;
#endif
}

iint CHecMPI::Send(void* buf, iint count, MPI_Datatype datatype, iint dest, iint tag, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Send(buf, count, datatype, dest, tag, comm);
#else
    return 0;
#endif
}
iint CHecMPI::Recv(void* buf, iint count, MPI_Datatype datatype, iint source, iint tag, MPI_Comm comm, MPI_Status *status)
{
#ifdef HAVE_MPI
    return MPI_Recv(buf, count, datatype, source, tag, comm, status);
#else
    return 0;
#endif
}

iint CHecMPI::Allgather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype,
                        void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
iint CHecMPI::Gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype,
                     void* recvbuf, iint recvcount, MPI_Datatype recvtype,
                     iint root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm);
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
iint CHecMPI::Scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype,
                      void* recvbuf, iint recvcnt, MPI_Datatype recvtype,
                      iint root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
#else
    recvbuf = sendbuf;
    return 0;
#endif
}
iint CHecMPI::Bcast(void* buf, iint cnt, MPI_Datatype type, iint root, MPI_Comm comm)
{
#ifdef HAVE_MPI
    return MPI_Bcast(buf, cnt, type, root, comm);
#else
    buf = buf;
    return 0;
#endif
}

MPI_Datatype CHecMPI::MPI_IINT()
{
    return mnMPI_IINT;
}
MPI_Datatype CHecMPI::MPI_UIINT()
{
    return mnMPI_UIINT;
}



