//
//
//  MW_Str.cpp
//  Structure Analysis "main" Routine
//
//  with pMW static lib
//
//          2009.7.29
//          2009.4.02
//          k.Takeda


#include "HEC_MW3.h"
using namespace pmw;

//#include "Logger.h"
//#include "BoundaryFace.h"


int main(int argc, char** argv)
{
//    string sInputFileName = "test_model_fukugou.msh";
//    string sOutputFileName= "test_model_fukugou.out";
    
//    string sInputFileName = "test_model_small_hexa_single.msh";
//    string sOutputFileName= "test_model_small_hexa_single.out";
//    string sInputFileName = "test_model_small_tetra_single.msh";
//    string sOutputFileName= "test_model_small_tetra_single.out";

//    string sInputFileName = "test_model_small_hexa.msh";
//    string sOutputFileName= "test_model_small_hexa.out";
    string sInputFileName = "test_model_small_tetra.msh";
    string sOutputFileName= "test_model_small_tetra.out";

//    string sInputFileName = "test_model_small_prism_single.msh";
//    string sOutputFileName= "test_model_small_prism_single.out";
//    string sInputFileName = "test_model_small_pyramid_single.msh";
//    string sOutputFileName= "test_model_small_pyramid_single.out";
    
//    string sInputFileName = "test_model_small_quad_single.msh";
//    string sOutputFileName= "test_model_small_quad_single.out";
//    string sInputFileName = "test_model_small_triangle_single.msh";
//    string sOutputFileName= "test_model_small_triangle_single.out";
//    string sInputFileName = "test_model_small_beam_single.msh";
//    string sOutputFileName= "test_model_small_beam_single.out";
    
    CMWMain *pMW = CMWMain::Instance();
    
    pMW->Initialize();
    //    pMW->InitGMG(5);     // Test Method <= OK.
    //    pMW->geneAssyModel();// Test Method <= OK.
    
    pMW->FileRead(sInputFileName);

    pMW->Refine();

    pMW->FileWrite(sOutputFileName);

    
    pMW->Finalize();

    return 0;
}
