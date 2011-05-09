//
//  EdgeFaceTree.cpp
//
//  Edge0, Edge1からFace局所インデックス番号を選択
//
//                          2009.06.30
//                          2009.06.30
//                          k.Takeda

#include "EdgeFaceTree.h"
#include "Logger.h"
using namespace pmw;

// construct
//
CEdgeFaceTree::CEdgeFaceTree()
{
}
// destruct
//
CEdgeFaceTree::~CEdgeFaceTree()
{
    cout << "~CEdgeFaceTree" << endl;
}

// method
// --
// 辺2本から面(Face)番号を出力
// --
// Hexa
//
uint& CEdgeFaceTree::getHexaFaceIndex(const uint& edge0, const uint& edge1)
{
    switch(edge0){
        case(0):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=0;}
            if(edge1==9){mFaceIndex=4;}
            if(edge1==8){mFaceIndex=4;}
            break;
        case(1):
            if(edge1==0) {mFaceIndex=0;}
            if(edge1==2) {mFaceIndex=0;}
            if(edge1==9) {mFaceIndex=2;}
            if(edge1==10){mFaceIndex=2;}
            break;
        case(2):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=0;}
            if(edge1==10){mFaceIndex=5;}
            if(edge1==11){mFaceIndex=5;}
            break;
        case(3):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==8){mFaceIndex=3;}
            if(edge1==11){mFaceIndex=3;}
            break;
        case(4):
            if(edge1==5){mFaceIndex=1;}
            if(edge1==7){mFaceIndex=1;}
            if(edge1==8){mFaceIndex=4;}
            if(edge1==9){mFaceIndex=4;}
            break;
        case(5):
            if(edge1==9){mFaceIndex=2;}
            if(edge1==10){mFaceIndex=2;}
            if(edge1==6){mFaceIndex=1;}
            if(edge1==4){mFaceIndex=1;}
            break;
        case(6):
            if(edge1==5){mFaceIndex=1;}
            if(edge1==7){mFaceIndex=1;}
            if(edge1==10){mFaceIndex=5;}
            if(edge1==11){mFaceIndex=5;}
            break;
        case(7):
            if(edge1==4){mFaceIndex=1;}
            if(edge1==6){mFaceIndex=1;}
            if(edge1==8){mFaceIndex=3;}
            if(edge1==11){mFaceIndex=3;}
            break;
        case(8):
            if(edge1==0){mFaceIndex=4;}
            if(edge1==4){mFaceIndex=4;}
            if(edge1==3){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            break;
        case(9):
            if(edge1==0){mFaceIndex=4;}
            if(edge1==4){mFaceIndex=4;}
            if(edge1==1){mFaceIndex=2;}
            if(edge1==5){mFaceIndex=2;}
            break;
        case(10):
            if(edge1==1){mFaceIndex=2;}
            if(edge1==5){mFaceIndex=2;}
            if(edge1==6){mFaceIndex=5;}
            if(edge1==2){mFaceIndex=5;}
            break;
        case(11):
            if(edge1==3){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            if(edge1==2){mFaceIndex=5;}
            if(edge1==6){mFaceIndex=5;}
            break;
        default:
            break;
    }

    return mFaceIndex;
}

// Tetra
// --
uint& CEdgeFaceTree::getTetraFaceIndex(const uint& edge0, const uint& edge1)
{
    switch(edge0){
        case(0):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=1;}
            if(edge1==4){mFaceIndex=1;}
            break;
        case(1):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==4){mFaceIndex=2;}
            if(edge1==5){mFaceIndex=2;}
            break;
        case(2):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==1){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=3;}
            if(edge1==5){mFaceIndex=3;}
            break;
        case(3):
            if(edge1==0){mFaceIndex=1;}
            if(edge1==4){mFaceIndex=1;}
            if(edge1==2){mFaceIndex=3;}
            if(edge1==5){mFaceIndex=3;}
            break;
        case(4):
            if(edge1==0){mFaceIndex=1;}
            if(edge1==3){mFaceIndex=1;}
            if(edge1==1){mFaceIndex=2;}
            if(edge1==5){mFaceIndex=2;}
            break;
        case(5):
            if(edge1==2){mFaceIndex=3;}
            if(edge1==3){mFaceIndex=3;}
            if(edge1==4){mFaceIndex=2;}
            if(edge1==1){mFaceIndex=2;}
            break;

        default:
            break;
    }

    return mFaceIndex;
}

// Prism
// --
uint& CEdgeFaceTree::getPrismFaceIndex(const uint& edge0, const uint& edge1)
{
    switch(edge0){
        case(0):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=2;}
            if(edge1==4){mFaceIndex=2;}
            break;
        case(1):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=4;}
            if(edge1==5){mFaceIndex=4;}
            break;
        case(2):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==1){mFaceIndex=0;}
            if(edge1==4){mFaceIndex=3;}
            if(edge1==5){mFaceIndex=3;}
            break;
        case(3):
            if(edge1==1){mFaceIndex=4;}
            if(edge1==8){mFaceIndex=4;}
            if(edge1==0){mFaceIndex=2;}
            if(edge1==6){mFaceIndex=2;}
            break;
        case(4):
            if(edge1==0){mFaceIndex=2;}
            if(edge1==6){mFaceIndex=2;}
            if(edge1==2){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            break;
        case(5):
            if(edge1==1){mFaceIndex=4;}
            if(edge1==8){mFaceIndex=4;}
            if(edge1==2){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            break;
        case(6):
            if(edge1==7){mFaceIndex=1;}
            if(edge1==8){mFaceIndex=1;}
            if(edge1==3){mFaceIndex=2;}
            if(edge1==4){mFaceIndex=2;}
            break;
        case(7):
            if(edge1==6){mFaceIndex=1;}
            if(edge1==8){mFaceIndex=1;}
            if(edge1==4){mFaceIndex=3;}
            if(edge1==5){mFaceIndex=3;}
            break;
        case(8):
            if(edge1==6){mFaceIndex=1;}
            if(edge1==7){mFaceIndex=1;}
            if(edge1==3){mFaceIndex=4;}
            if(edge1==5){mFaceIndex=4;}
            break;

        default:
            break;
    }

    return mFaceIndex;
}

// Pyramid
// --
uint& CEdgeFaceTree::getPyramidFaceIndex(const uint& edge0, const uint& edge1)
{
    switch(edge0){
        case(0):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=0;}
            if(edge1==4){mFaceIndex=4;}
            if(edge1==7){mFaceIndex=4;}
            break;
        case(1):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==4){mFaceIndex=1;}
            if(edge1==5){mFaceIndex=1;}
            break;
        case(2):
            if(edge1==1){mFaceIndex=0;}
            if(edge1==3){mFaceIndex=0;}
            if(edge1==5){mFaceIndex=1;}
            if(edge1==6){mFaceIndex=1;}
            break;
        case(3):
            if(edge1==0){mFaceIndex=0;}
            if(edge1==2){mFaceIndex=0;}
            if(edge1==6){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            break;
        case(4):
            if(edge1==0){mFaceIndex=4;}
            if(edge1==7){mFaceIndex=4;}
            if(edge1==1){mFaceIndex=1;}
            if(edge1==5){mFaceIndex=1;}
            break;
        case(5):
            if(edge1==1){mFaceIndex=1;}
            if(edge1==4){mFaceIndex=1;}
            if(edge1==2){mFaceIndex=2;}
            if(edge1==6){mFaceIndex=2;}
            break;
        case(6):
            if(edge1==2){mFaceIndex=2;}
            if(edge1==5){mFaceIndex=2;}
            if(edge1==3){mFaceIndex=3;}
            if(edge1==7){mFaceIndex=3;}
            break;
        case(7):
            if(edge1==0){mFaceIndex=4;}
            if(edge1==4){mFaceIndex=4;}
            if(edge1==3){mFaceIndex=3;}
            if(edge1==6){mFaceIndex=3;}
            break;

        default:
            break;
    }

    return mFaceIndex;
}

// Quad
// --
uint& CEdgeFaceTree::getQuadFaceIndex(const uint& edge0, const uint& edge1)
{
    if(edge0 > 3 || edge1 > 3){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"edge number over flow @CEdgeFaceTree::getQuadFaceIndex");
    }

    mFaceIndex=0;

    return mFaceIndex;
}

// Triangle
// --
uint& CEdgeFaceTree::getTriangleFaceIndex(const uint& edge0, const uint& edge1)
{
    if(edge0 > 2 || edge1 > 2){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"edge number over flow @CEdgeFaceTree::getTriangleFaceIndex");
    }
    
    mFaceIndex=0;

    return mFaceIndex;
}










