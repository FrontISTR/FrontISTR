/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/HEC_MPI.h
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
#define MPI_CHAR 0
#define MPI_WCHAR 0
#define MPI_SHORT 0
#define MPI_INT 0
#define MPI_LONG 0
#define MPI_SIGNED_CHAR 0
#define MPI_UNSIGNED_CHAR 0
#define MPI_UNSIGNED_SHORT 0
#define MPI_UNSIGNED 0        
#define MPI_UNSIGEND_LONG 0
#define MPI_FLOAT 0
#define MPI_DOUBLE 0
#define MPI_SUM 0
#define MPI_PROD 0
#define MPI_MAX 0
#define MPI_MIN 0
#define MPI_COMM_WORLD 0
#endif
#include <iostream>
using namespace std;
#include "TypeDef.h"
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
    int nNumOfProcess;
public:
    void Initialize(int argc, char **argv);
    void Finalize();
    int& getRank(){ return myRank;}
    int& getNumOfProcess(){ return nNumOfProcess;}
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
    int Bcast(void* buf, int cnt, MPI_Datatype type, int root, MPI_Comm comm);
};
#endif	/* _HEC_MPI_H */
}
