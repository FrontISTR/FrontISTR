//
//  HEC_MW.cpp
//
//  mpichのラッパー
//
//                  2009.08.24
//                  2009.08.24
//                  k.Takeda
#include "HEC_MPI.h"
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

// 
//
void CHecMPI::Initialize(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //debug
    std::cout << " myRank => " << myRank << ",number of process => " << numOfProcess << std::endl;
}

void CHecMPI::Finalize()
{
    MPI_Finalize();
    
}






