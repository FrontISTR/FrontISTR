/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/EdgeFaceTree.cpp
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
#include "EdgeFaceTree.h"
#include "Logger.h"
using namespace pmw;
CEdgeFaceTree::CEdgeFaceTree()
{
    uiint iface,iedge;
    uiint hexedge[6][4]= {
        {0,1,2,3},
        {4,7,6,5},
        {1,9,5,10},
        {11,7,8,3},
        {0,8.4,9},
        {2,10,6,11}
    };
    for(iface=0; iface< 6; iface++) {
        for(iedge=0; iedge< 4; iedge++) {
            mHexaFaceConnEdge[iface][iedge]= hexedge[iface][iedge];
        };
    };
    uiint tetedge[4][3]= {
        {0,1,2},
        {0,3,4},
        {4,5,1},
        {2,5,3}
    };
    for(iface=0; iface< 4; iface++) {
        for(iedge=0; iedge< 3; iedge++) {
            mTetraFaceConnEdge[iface][iedge]= tetedge[iface][iedge];
        };
    };
    uiint priedge[5][4]= {
        {0,2,1,0},
        {8,7,6,8},
        {0,3,6,4},
        {4,7,5,2},
        {1,5,8,3}
    };
    for(iface=0; iface< 5; iface++) {
        for(iedge=0; iedge< 4; iedge++) {
            mPrismFaceConnEdge[iface][iedge]= priedge[iface][iedge];
        };
    };
    uiint pyraedge[5][4]= {
        {0,1,2,3},
        {1,4,5,1},
        {2,5,6,2},
        {6,1,3,6},
        {0,7,4,0}
    };
    for(iface=0; iface< 5; iface++) {
        for(iedge=0; iedge< 4; iedge++) {
            mPyramidFaceConnEdge[iface][iedge]= pyraedge[iface][iedge];
        };
    };
    uiint quadedge[1][4]= {
        {0,1,2,3}
    };
    for(iedge=0; iedge< 4; iedge++) {
        mQuadFaceConnEdge[0][iedge]= quadedge[0][iedge];
    }
    uiint triedge[1][3]= {
        {0,1,2}
    };
    for(iedge=0; iedge< 3; iedge++) {
        mTriangleFaceConnEdge[0][iedge]= triedge[0][iedge];
    }
}
CEdgeFaceTree::~CEdgeFaceTree()
{
    cout << "~CEdgeFaceTree" << endl;
}
uiint& CEdgeFaceTree::getHexaFaceIndex(const uiint& edge0, const uiint& edge1)
{
    switch(edge0) {
    case(0):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=0;
        }
        if(edge1==9) {
            mFaceIndex=4;
        }
        if(edge1==8) {
            mFaceIndex=4;
        }
        break;
    case(1):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==9) {
            mFaceIndex=2;
        }
        if(edge1==10) {
            mFaceIndex=2;
        }
        break;
    case(2):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=0;
        }
        if(edge1==10) {
            mFaceIndex=5;
        }
        if(edge1==11) {
            mFaceIndex=5;
        }
        break;
    case(3):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==8) {
            mFaceIndex=3;
        }
        if(edge1==11) {
            mFaceIndex=3;
        }
        break;
    case(4):
        if(edge1==5) {
            mFaceIndex=1;
        }
        if(edge1==7) {
            mFaceIndex=1;
        }
        if(edge1==8) {
            mFaceIndex=4;
        }
        if(edge1==9) {
            mFaceIndex=4;
        }
        break;
    case(5):
        if(edge1==9) {
            mFaceIndex=2;
        }
        if(edge1==10) {
            mFaceIndex=2;
        }
        if(edge1==6) {
            mFaceIndex=1;
        }
        if(edge1==4) {
            mFaceIndex=1;
        }
        break;
    case(6):
        if(edge1==5) {
            mFaceIndex=1;
        }
        if(edge1==7) {
            mFaceIndex=1;
        }
        if(edge1==10) {
            mFaceIndex=5;
        }
        if(edge1==11) {
            mFaceIndex=5;
        }
        break;
    case(7):
        if(edge1==4) {
            mFaceIndex=1;
        }
        if(edge1==6) {
            mFaceIndex=1;
        }
        if(edge1==8) {
            mFaceIndex=3;
        }
        if(edge1==11) {
            mFaceIndex=3;
        }
        break;
    case(8):
        if(edge1==0) {
            mFaceIndex=4;
        }
        if(edge1==4) {
            mFaceIndex=4;
        }
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        break;
    case(9):
        if(edge1==0) {
            mFaceIndex=4;
        }
        if(edge1==4) {
            mFaceIndex=4;
        }
        if(edge1==1) {
            mFaceIndex=2;
        }
        if(edge1==5) {
            mFaceIndex=2;
        }
        break;
    case(10):
        if(edge1==1) {
            mFaceIndex=2;
        }
        if(edge1==5) {
            mFaceIndex=2;
        }
        if(edge1==6) {
            mFaceIndex=5;
        }
        if(edge1==2) {
            mFaceIndex=5;
        }
        break;
    case(11):
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        if(edge1==2) {
            mFaceIndex=5;
        }
        if(edge1==6) {
            mFaceIndex=5;
        }
        break;
    default:
        break;
    }
    return mFaceIndex;
}
uiint& CEdgeFaceTree::getTetraFaceIndex(const uiint& edge0, const uiint& edge1)
{
    switch(edge0) {
    case(0):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=1;
        }
        if(edge1==4) {
            mFaceIndex=1;
        }
        break;
    case(1):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==4) {
            mFaceIndex=2;
        }
        if(edge1==5) {
            mFaceIndex=2;
        }
        break;
    case(2):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==5) {
            mFaceIndex=3;
        }
        break;
    case(3):
        if(edge1==0) {
            mFaceIndex=1;
        }
        if(edge1==4) {
            mFaceIndex=1;
        }
        if(edge1==2) {
            mFaceIndex=3;
        }
        if(edge1==5) {
            mFaceIndex=3;
        }
        break;
    case(4):
        if(edge1==0) {
            mFaceIndex=1;
        }
        if(edge1==3) {
            mFaceIndex=1;
        }
        if(edge1==1) {
            mFaceIndex=2;
        }
        if(edge1==5) {
            mFaceIndex=2;
        }
        break;
    case(5):
        if(edge1==2) {
            mFaceIndex=3;
        }
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==4) {
            mFaceIndex=2;
        }
        if(edge1==1) {
            mFaceIndex=2;
        }
        break;
    default:
        break;
    }
    return mFaceIndex;
}
uiint& CEdgeFaceTree::getPrismFaceIndex(const uiint& edge0, const uiint& edge1)
{
    switch(edge0) {
    case(0):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=2;
        }
        if(edge1==4) {
            mFaceIndex=2;
        }
        break;
    case(1):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=4;
        }
        if(edge1==5) {
            mFaceIndex=4;
        }
        break;
    case(2):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==4) {
            mFaceIndex=3;
        }
        if(edge1==5) {
            mFaceIndex=3;
        }
        break;
    case(3):
        if(edge1==1) {
            mFaceIndex=4;
        }
        if(edge1==8) {
            mFaceIndex=4;
        }
        if(edge1==0) {
            mFaceIndex=2;
        }
        if(edge1==6) {
            mFaceIndex=2;
        }
        break;
    case(4):
        if(edge1==0) {
            mFaceIndex=2;
        }
        if(edge1==6) {
            mFaceIndex=2;
        }
        if(edge1==2) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        break;
    case(5):
        if(edge1==1) {
            mFaceIndex=4;
        }
        if(edge1==8) {
            mFaceIndex=4;
        }
        if(edge1==2) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        break;
    case(6):
        if(edge1==7) {
            mFaceIndex=1;
        }
        if(edge1==8) {
            mFaceIndex=1;
        }
        if(edge1==3) {
            mFaceIndex=2;
        }
        if(edge1==4) {
            mFaceIndex=2;
        }
        break;
    case(7):
        if(edge1==6) {
            mFaceIndex=1;
        }
        if(edge1==8) {
            mFaceIndex=1;
        }
        if(edge1==4) {
            mFaceIndex=3;
        }
        if(edge1==5) {
            mFaceIndex=3;
        }
        break;
    case(8):
        if(edge1==6) {
            mFaceIndex=1;
        }
        if(edge1==7) {
            mFaceIndex=1;
        }
        if(edge1==3) {
            mFaceIndex=4;
        }
        if(edge1==5) {
            mFaceIndex=4;
        }
        break;
    default:
        break;
    }
    return mFaceIndex;
}
uiint& CEdgeFaceTree::getPyramidFaceIndex(const uiint& edge0, const uiint& edge1)
{
    switch(edge0) {
    case(0):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=0;
        }
        if(edge1==4) {
            mFaceIndex=4;
        }
        if(edge1==7) {
            mFaceIndex=4;
        }
        break;
    case(1):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==4) {
            mFaceIndex=1;
        }
        if(edge1==5) {
            mFaceIndex=1;
        }
        break;
    case(2):
        if(edge1==1) {
            mFaceIndex=0;
        }
        if(edge1==3) {
            mFaceIndex=0;
        }
        if(edge1==5) {
            mFaceIndex=1;
        }
        if(edge1==6) {
            mFaceIndex=1;
        }
        break;
    case(3):
        if(edge1==0) {
            mFaceIndex=0;
        }
        if(edge1==2) {
            mFaceIndex=0;
        }
        if(edge1==6) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        break;
    case(4):
        if(edge1==0) {
            mFaceIndex=4;
        }
        if(edge1==7) {
            mFaceIndex=4;
        }
        if(edge1==1) {
            mFaceIndex=1;
        }
        if(edge1==5) {
            mFaceIndex=1;
        }
        break;
    case(5):
        if(edge1==1) {
            mFaceIndex=1;
        }
        if(edge1==4) {
            mFaceIndex=1;
        }
        if(edge1==2) {
            mFaceIndex=2;
        }
        if(edge1==6) {
            mFaceIndex=2;
        }
        break;
    case(6):
        if(edge1==2) {
            mFaceIndex=2;
        }
        if(edge1==5) {
            mFaceIndex=2;
        }
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==7) {
            mFaceIndex=3;
        }
        break;
    case(7):
        if(edge1==0) {
            mFaceIndex=4;
        }
        if(edge1==4) {
            mFaceIndex=4;
        }
        if(edge1==3) {
            mFaceIndex=3;
        }
        if(edge1==6) {
            mFaceIndex=3;
        }
        break;
    default:
        break;
    }
    return mFaceIndex;
}
uiint& CEdgeFaceTree::getQuadFaceIndex(const uiint& edge0, const uiint& edge1)
{
    if(edge0 > 3 || edge1 > 3) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"edge number over flow @CEdgeFaceTree::getQuadFaceIndex");
    }
    mFaceIndex=0;
    return mFaceIndex;
}
uiint& CEdgeFaceTree::getTriangleFaceIndex(const uiint& edge0, const uiint& edge1)
{
    if(edge0 > 2 || edge1 > 2) {
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"edge number over flow @CEdgeFaceTree::getTriangleFaceIndex");
    }
    mFaceIndex=0;
    return mFaceIndex;
}
