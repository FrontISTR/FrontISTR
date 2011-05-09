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

//#include "/usr/lib/mpich/include/mpi.h"
#include "mpi/mpi.h"

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
};
#endif	/* _HEC_MPI_H */
}
