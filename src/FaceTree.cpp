/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FaceTree.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FaceTree.h"
#include "Logger.h"
using namespace pmw;
CFaceTree::CFaceTree()
{
    uint i,ii;    
    uint hface[6][4]={
    {0,1,2,3},
    {4,7,6,5},
    {1,5,6,2},
    {7,4,0,3},
    {0,4,5,1},
    {2,6,7,3}
    };
    for(i=0; i< 6; i++){
    for(ii=0; ii< 4; ii++){
        mHexaFaceIndex[i][ii]= hface[i][ii];
    };
    };
    uint tface[4][3]={
    {0,1,2},
    {0,3,1},
    {1,3,2},
    {0,2,3}
    };
    for(i=0; i< 4; i++){
    for(ii=0; ii< 3; ii++){
       mTetraFaceIndex[i][ii]= tface[i][ii];
    };
    };
    uint pface[5][4]={
    {0,1,2,0},
    {3,5,4,3},
    {0,3,4,1},
    {1,2,5,4},
    {0,2,5,3} 
    };
    for(i=0; i< 5; i++){
    for(ii=0; ii< 4; ii++){
       mPrismFaceIndex[i][ii]= pface[i][ii];
    };
    };
    uint pyface[5][4]={
    {0,1,2,3},
    {1,2,4,1},
    {2,4,3,2},
    {3,4,0,3},
    {0,4,1,0} 
    };
    for(i=0; i< 5; i++){
    for(ii=0; ii< 4; ii++){
       mPyramidFaceIndex[i][ii]= pyface[i][ii];
    };
    };
    uint quad[4]={0,1,2,3};
    for(i=0; i< 4; i++){
        mQuadFaceIndex[i]= quad[i];
    };
    uint tri[3]={0,1,2};
    for(i=0; i< 3; i++){
        mTriangleFaceIndex[i]= tri[i];
    };
    uint hexconn[8][3]={
        {0,4,3},
        {0,4,2},
        {0,2,5},
        {0,3,5},
        {1,4,3},
        {1,4,2},
        {1,2,5},
        {1,3,5}
    };
    for(i=0; i< 8; i++){
        for(ii=0; ii< 3; ii++){
            mHexaConnFace[i][ii]= hexconn[i][ii];
        }
    }
    uint tetconn[4][3]={
        {0,3,1},
        {0,1,2},
        {0,2,3},
        {1,2,3}
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 3; ii++){
            mTetraConnFace[i][ii]= tetconn[i][ii];
        }
    }
    uint priconn[6][3]={
        {0,4,2},
        {0,2,3},
        {0,4,3},
        {1,4,2},
        {1,2,3},
        {1,3,4}
    };
    for(i=0; i< 6; i++){
        for(ii=0; ii< 3; ii++){
            mPrismConnFace[i][ii]= priconn[i][ii];
        }
    }
    uint pyconn[5][4]={
        {0,4,3,0},
        {0,4,1,0},
        {0,1,2,0},
        {0,2,3,0},
        {1,2,3,4}
    };
    for(i=0; i< 5; i++){
        for(ii=0; ii< 4; ii++){
            mPyramidConnFace[i][ii]= pyconn[i][ii];
        }
    }
    uint quconn[4][1]={
        {0},{0},{0},{0}
    };
    for(i=0; i< 4; i++){
        mQuadConnFace[i][0]= quconn[i][0];
    }
    uint triconn[3][1]={
        {0},{0},{0}
    };
    for(i=0; i< 3; i++){
        mTriangleConnFace[i][0]= triconn[i][0];
    }
}
CFaceTree::~CFaceTree()
{
}
uint& CFaceTree::getHexaFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;
    uint ia,ib;
    for(faceIndex=0; faceIndex< 6; faceIndex++){
        hitCount=0;
        for(ia=0; ia< 4; ia++){
        for(ib=0; ib< 4; ib++){
            if(mHexaFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]){ ++hitCount; break;}
        };
        };
        if(hitCount==4){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
uint& CFaceTree::getHexaFaceIndex2(const vuint& vLocalNodeIndex)
{
    switch(vLocalNodeIndex[0]){
    case(0):
        switch(vLocalNodeIndex[2]){
            case(2): mFaceIndex=0; break;
            case(5): mFaceIndex=4; break;
            case(7): mFaceIndex=3; break;           
        }
        break;
    case(1):
        switch(vLocalNodeIndex[2]){
            case(3): mFaceIndex=0; break;
            case(6): mFaceIndex=2; break;
            case(4): mFaceIndex=4; break;
        }
        break;
    case(2):
        switch(vLocalNodeIndex[2]){
            case(5): mFaceIndex=2; break;
            case(0): mFaceIndex=0; break;
            case(7): mFaceIndex=5; break;
        }
        break;
    case(3):
        switch(vLocalNodeIndex[2]){
            case(1): mFaceIndex=0; break;
            case(6): mFaceIndex=5; break;
            case(4): mFaceIndex=3; break;
        }
        break;
    case(4):
        switch(vLocalNodeIndex[2]){
            case(1): mFaceIndex=4; break;
            case(3): mFaceIndex=3; break;
            case(6): mFaceIndex=1; break;
        }
        break;
    case(5):
        switch(vLocalNodeIndex[2]){
            case(2): mFaceIndex=2; break;
            case(0): mFaceIndex=4; break;
            case(7): mFaceIndex=1; break;
        }
        break;
    case(6):
        switch(vLocalNodeIndex[2]){
            case(1): mFaceIndex=2; break;
            case(3): mFaceIndex=5; break;
            case(4): mFaceIndex=1; break;
        }
        break;
    case(7):
        switch(vLocalNodeIndex[2]){
            case(2): mFaceIndex=5; break;
            case(5): mFaceIndex=1; break;
            case(0): mFaceIndex=3; break;
        }
        break;
    }
    return mFaceIndex;
}
uint& CFaceTree::getTetraFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;
    uint ia,ib;
    for(faceIndex=0; faceIndex< 4; faceIndex++){
        hitCount=0;
        for(ia=0; ia< 3; ia++){
        for(ib=0; ib< 3; ib++){
            if(mTetraFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]) ++hitCount;
        };
        };
        if(hitCount==3){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
uint& CFaceTree::getTetraFaceIndex2(const vuint& vLocalNodeIndex)
{
    switch(vLocalNodeIndex[0]){
        case(0):
            switch(vLocalNodeIndex[1]){
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(2): mFaceIndex=0; break;
                    case(3): mFaceIndex=1; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(1): mFaceIndex=1; break;
                    case(2): mFaceIndex=3; break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(1): mFaceIndex=0; break;
                    case(3): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(1):
            switch(vLocalNodeIndex[1]){
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(2): mFaceIndex=0; break;
                    case(3): mFaceIndex=1; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=1; break;
                    case(2): mFaceIndex=2; break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=0; break;
                    case(3): mFaceIndex=2; break;
                    }
                    break;
            }
            break;
        case(2):
            switch(vLocalNodeIndex[1]){
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=3; break;
                    case(1): mFaceIndex=2; break;
                    }
                    break;
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=0; break;
                    case(3): mFaceIndex=2; break;
                    }
                    break;
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=3; break;
                    case(1): mFaceIndex=0; break;
                    }
                    break;
            }
            break;
        case(3):
            switch(vLocalNodeIndex[1]){
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=3; break;
                    case(1): mFaceIndex=2; break;
                    }
                    break;
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(2): mFaceIndex=3; break;
                    case(1): mFaceIndex=1; break;
                    }
                    break;
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=1; break;
                    case(2): mFaceIndex=2; break;
                    }
                    break;
            }
            break;
    }
    return mFaceIndex;
}
uint& CFaceTree::getPrismFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;
    uint ia,ib;
    for(faceIndex=0; faceIndex< 5; faceIndex++){
        hitCount=0;
        for(ia=0; ia< 3; ia++){
        for(ib=0; ib< 3; ib++){
            if(mPrismFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]) ++hitCount;
        };
        };
        if(hitCount==3){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
uint& CFaceTree::getPrismFaceIndex2(const vuint& vLocalNodeIndex)
{
    switch(vLocalNodeIndex[0]){
        case(0):
            switch(vLocalNodeIndex[1]){
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(4): mFaceIndex=2; break;
                    case(2): mFaceIndex=0; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(4): mFaceIndex=2; break;
                    case(5): mFaceIndex=4; break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(1): mFaceIndex=0; break;
                    case(5): mFaceIndex=4; break;
                    }
                    break;
            }
            break;
        case(1):
            switch(vLocalNodeIndex[1]){
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=2; break;
                    case(2): mFaceIndex=0; break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=0; break;
                    case(5): mFaceIndex=3; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=2; break;
                    case(5): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(2):
            switch(vLocalNodeIndex[1]){
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(1): mFaceIndex=0; break;
                    case(3): mFaceIndex=4; break;
                    }
                    break;
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=0; break;
                    case(4): mFaceIndex=3; break;
                    }
                    break;
                case(5):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=4; break;
                    case(4): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(3):
            switch(vLocalNodeIndex[1]){
                case(5):
                    switch(vLocalNodeIndex[2]){
                    case(2): mFaceIndex=4; break;
                    case(4): mFaceIndex=1; break;
                    }
                    break;
                case(0):
                    switch(vLocalNodeIndex[2]){
                    case(1): mFaceIndex=2; break;
                    case(2): mFaceIndex=4; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                    case(5): mFaceIndex=1; break;
                    case(1): mFaceIndex=2; break;
                    }
                    break;
            }
            break;
        case(4):
            switch(vLocalNodeIndex[1]){
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(5): mFaceIndex=1; break;
                    case(0): mFaceIndex=2; break;
                    }
                    break;
                case(1):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=2; break;
                    case(2): mFaceIndex=3; break;
                    }
                    break;
                case(5):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=1; break;
                    case(2): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(5):
            switch(vLocalNodeIndex[1]){
                case(2):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=4; break;
                    case(1): mFaceIndex=3; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                    case(0): mFaceIndex=4; break;
                    case(4): mFaceIndex=1; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                    case(3): mFaceIndex=1; break;
                    case(1): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
    }
    return mFaceIndex;
}
uint& CFaceTree::getPyramidFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;
    uint ia,ib;
    for(faceIndex=0; faceIndex< 5; faceIndex++){
        hitCount=0;
        for(ia=0; ia< 3; ia++){
        for(ib=0; ib< 3; ib++){
            if(mPyramidFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]) ++hitCount;
        };
        };
        if(hitCount==3){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
uint& CFaceTree::getPyramidFaceIndex2(const vuint& vLocalNodeIndex)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    switch(vLocalNodeIndex[0]){
        case(0):
            switch(vLocalNodeIndex[1]){
                case(1):
                    switch(vLocalNodeIndex[2]){
                        case(2): mFaceIndex=0; break;
                        case(4): mFaceIndex=4; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                        case(2): mFaceIndex=0; break;
                        case(4): mFaceIndex=3; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=4; break;
                        case(2): mFaceIndex=5;pLogger->Info(Utility::LoggerMode::Error, "CFaceTree::getPyramidFaceIndex2 case:0-4-2"); break;
                        case(3): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(1):
            switch(vLocalNodeIndex[1]){
                case(0):
                    switch(vLocalNodeIndex[2]){
                        case(3): mFaceIndex=0; break;
                        case(4): mFaceIndex=4; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                        case(0): mFaceIndex=4; break;
                        case(2): mFaceIndex=1; break;
                        case(3): mFaceIndex=5;pLogger->Info(Utility::LoggerMode::Error, "CFaceTree::getPyramidFaceIndex2 case:1-4-3"); break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                        case(3): mFaceIndex=0; break;
                        case(4): mFaceIndex=1; break;
                    }
                    break;
            }
            break;
        case(2):
            switch(vLocalNodeIndex[1]){
                case(1):
                    switch(vLocalNodeIndex[2]){
                        case(4): mFaceIndex=1; break;
                        case(0): mFaceIndex=0; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                        case(0): mFaceIndex=0; break;
                        case(4): mFaceIndex=2; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=1; break;
                        case(0): mFaceIndex=5;pLogger->Info(Utility::LoggerMode::Error, "CFaceTree::getPyramidFaceIndex2 case:2-4-0"); break;
                        case(3): mFaceIndex=2; break;
                    }
                    break;
            }
            break;
        case(3):
            switch(vLocalNodeIndex[1]){
                case(2):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=0; break;
                        case(4): mFaceIndex=2; break;
                    }
                    break;
                case(4):
                    switch(vLocalNodeIndex[2]){
                        case(2): mFaceIndex=2; break;
                        case(1): mFaceIndex=5;pLogger->Info(Utility::LoggerMode::Error, "CFaceTree::getPyramidFaceIndex2 case:3-4-1"); break;
                        case(0): mFaceIndex=3; break;
                    }
                    break;
                case(0):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=0; break;
                        case(4): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
        case(4):
            switch(vLocalNodeIndex[1]){
                case(0):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=4; break;
                        case(3): mFaceIndex=3; break;
                    }
                    break;
                case(1):
                    switch(vLocalNodeIndex[2]){
                        case(0): mFaceIndex=4; break;
                        case(2): mFaceIndex=1; break;
                    }
                    break;
                case(2):
                    switch(vLocalNodeIndex[2]){
                        case(1): mFaceIndex=1; break;
                        case(3): mFaceIndex=2; break;
                    }
                    break;
                case(3):
                    switch(vLocalNodeIndex[2]){
                        case(2): mFaceIndex=2; break;
                        case(0): mFaceIndex=3; break;
                    }
                    break;
            }
            break;
    }
    return mFaceIndex;
}
uint& CFaceTree::getQuadFaceIndex(const vuint& vLocalNodeIndex)
{
    uint ia,ib;
    uint hitCount=0;
    uint numOfLocalNode = vLocalNodeIndex.size();
    for(ia=0; ia< 4; ia++){
    for(ib=0; ib< numOfLocalNode; ib++){
        if(mQuadFaceIndex[ia]==vLocalNodeIndex[ib]) ++hitCount;
    };
    };
    if(hitCount >= numOfLocalNode){
        mFaceIndex=0;
    }else{
        mFaceIndex=1;
    }
    return mFaceIndex;
}
uint& CFaceTree::getTriangleFaceIndex(const vuint& vLocalNodeIndex)
{
    uint ia,ib;
    uint hitCount=0;
    uint numOfLocalNode = vLocalNodeIndex.size();
    for(ia=0; ia< 3; ia++){
    for(ib=0; ib< numOfLocalNode; ib++){
        if(mTriangleFaceIndex[ia]==vLocalNodeIndex[ib]) ++hitCount;
    };
    };
    if(hitCount==numOfLocalNode){
        mFaceIndex=0;
    }else{
        mFaceIndex=1;
    }
    return mFaceIndex;
}
uint* CFaceTree::getLocalNodeHexaFace(const uint& faceIndex)
{
    if(faceIndex >= 6){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodeHexaFace",faceIndex);
    }
    return mHexaFaceIndex[faceIndex];
}
uint* CFaceTree::getLocalNodeTetraFace(const uint& faceIndex)
{
    if(faceIndex >= 4){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodeTetraFace",faceIndex);
    }
    return mTetraFaceIndex[faceIndex];
}
uint* CFaceTree::getLocalNodePrismFace(const uint& faceIndex)
{
    if(faceIndex >= 5){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodePrismFace",faceIndex);
    }
    return mPrismFaceIndex[faceIndex];
}
uint* CFaceTree::getLocalNodePyramidFace(const uint& faceIndex)
{
    if(faceIndex >= 5){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodePyramidFace",faceIndex);
    }
    return mPyramidFaceIndex[faceIndex];
}
uint* CFaceTree::getLocalNodeQuadFace()
{
    return mQuadFaceIndex;
}
uint* CFaceTree::getLocalNodeTriangleFace()
{
    return mTriangleFaceIndex;
}
