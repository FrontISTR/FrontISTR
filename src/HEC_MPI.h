/* 
 * File:   HEC_MPI.h
 * Author: ktakeda
 *
 * mpichラッパー
 *
 * Created on 2009/08/24, 14:13
 */
#include <iostream>
using namespace std;

#ifdef HAVE_MPI
#include "mpi.h"
#else
typedef int MPI_Comm;
typedef int MPI_Group;
typedef int MPI_Request;
typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Fint;

// MPI_Datatype
//
#define MPI_CHAR 0
#define MPI_WCHAR 0
#define MPI_SHORT 0
#define MPI_INT 0
#define MPI_LONG 0
#define MPI_SIGNED_CHAR 0
#define MPI_UNSIGNED_CHAR 0
#define MPI_UNSIGNED_SHORT 0
#define MPI_UNSIGNED 0        //unsigned int
#define MPI_UNSIGEND_LONG 0
#define MPI_FLOAT 0
#define MPI_DOUBLE 0
#define MPI_SUM 0
#define MPI_PROD 0
#define MPI_MAX 0
#define MPI_MIN 0
#define MPI_COMM_WORLD 0
#endif


namespace pmw{
#ifndef _HEC_MPI_H
#define	_HEC_MPI_H
class CHecMPI{
public:
    static CHecMPI* Instance(){
        static CHecMPI moMPI;
        return &moMPI;
    }
private:
    CHecMPI();
public:
    virtual ~CHecMPI();



private:
    int myRank;
    int numOfProcess;


public:
    void Initialize(int argc, char **argv);//引数:argc,argvは,MPIの引数
    void Finalize();

    int& getRank(){ return myRank;}//自分のプロセス-ランクを取得
    int& getNumOfProcess(){ return numOfProcess;}

    // MPI ラップ関数
    //
    int Isend(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
    int Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);

    int Wait(MPI_Request *request, MPI_Status *status);
    int Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
    int Barrier(MPI_Comm comm);
    int Abort(MPI_Comm comm, int errorcode);
    int Comm_rank(MPI_Comm comm, int *rank);

    int Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    int Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

    int Allgather(void* sendbuf, int sendcnt, MPI_Datatype sendtype,
                  void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    int Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm);
    int Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype,
                void* recvbuf, int recvcnt, MPI_Datatype recvtype,
                int root, MPI_Comm comm);

};
#endif	/* _HEC_MPI_H */
}










