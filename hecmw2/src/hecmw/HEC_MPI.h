/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/HEC_MPI.h
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
#ifndef _78eb95b5_2320_402e_a820_612db8db0f07
#define	_78eb95b5_2320_402e_a820_612db8db0f07

#include "CommonStd.h"
#include <iostream>
#include <sstream> //stringstream
using namespace std;
#include "TypeDef.h"

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


namespace pmw
{
class CHecMPI
{
public:
    static CHecMPI* Instance() {
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

    MPI_Datatype mnMPI_UIINT;
    MPI_Datatype mnMPI_IINT;

    //--
    // 派生データ型(IINT,UIINT) Allreduce
    //--
    iint Allreduce_sum(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm);
    iint Allreduce_max(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm);
    iint Allreduce_min(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm);
    iint Allreduce_prod(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Comm comm);

public:
    void Initialize(int argc, char **argv);
    void Finalize();
    int& getRank() {
        return myRank;
    }
    void Comm_rank(MPI_Comm comm, int *rank);
    int& getNumOfProcess() {
        return nNumOfProcess;
    }

    iint Isend(void* buf, iint count, MPI_Datatype datatype, iint dest, iint tag, MPI_Comm comm, MPI_Request *request);
    iint Irecv(void* buf, iint count, MPI_Datatype datatype, iint source, iint tag, MPI_Comm comm, MPI_Request *request);
    iint Wait(MPI_Request *request, MPI_Status *status);
    iint Waitall(iint nNumOfReq, MPI_Request *request, MPI_Status *status);

    iint Allreduce(void* sendbuf, void* recvbuf, iint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
    iint Barrier(MPI_Comm comm);
    iint Abort(MPI_Comm comm, iint errorcode);

    iint Send(void* buf, iint count, MPI_Datatype datatype, iint dest, iint tag, MPI_Comm comm);
    iint Recv(void* buf, iint count, MPI_Datatype datatype, iint source, iint tag, MPI_Comm comm, MPI_Status *status);

    iint Allgather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype,
                   void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    iint Gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype,
                void* recvbuf, iint recvcount, MPI_Datatype recvtype,
                iint root, MPI_Comm comm);
    iint Scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype,
                 void* recvbuf, iint recvcnt, MPI_Datatype recvtype,
                 iint root, MPI_Comm comm);
    iint Bcast(void* buf, iint cnt, MPI_Datatype type, iint root, MPI_Comm comm);

    MPI_Datatype MPI_UIINT();
    MPI_Datatype MPI_IINT();
};
}
#endif	 /* _78eb95b5_2320_402e_a820_612db8db0f07 */ /* _HEC_MPI_H */

