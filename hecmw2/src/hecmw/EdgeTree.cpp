/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/EdgeTree.cpp
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
#include "EdgeTree.h"
#include "ElementType.h"
#include "Logger.h"
#include "HEC_MPI.h"
using namespace pmw;
#include <iostream>
CEdgeTree::CEdgeTree()
{
    mpLogger = Utility::CLogger::Instance();
    uiint i,ii;
    uiint hexedge[12][2]= {
        {0,1},
        {1,2},
        {2,3},
        {3,0},
        {4,5},
        {5,6},
        {6,7},
        {7,4},
        {0,4},
        {1,5},
        {2,6},
        {3,7}
    };
    for(i=0; i< 12; i++) {
        for(ii=0; ii< 2; ii++) {
            mHexaEdgeIndex[i][ii]= hexedge[i][ii];
        };
    };
    uiint tetedge[6][2]= {
        {0,1},
        {1,2},
        {2,0},
        {0,3},
        {1,3},
        {2,3}
    };
    for(i=0; i< 6; i++) {
        for(ii=0; ii< 2; ii++) {
            mTetraEdgeIndex[i][ii]= tetedge[i][ii];
        };
    };
    uiint priedge[9][2]= {
        {0,1},
        {0,2},
        {1,2},
        {0,3},
        {1,4},
        {2,5},
        {3,4},
        {4,5},
        {5,3}
    };
    for(i=0; i< 9; i++) {
        for(ii=0; ii< 2; ii++) {
            mPrismEdgeIndex[i][ii]= priedge[i][ii];
        };
    };
    uiint pyedge[8][2]= {
        {0,1},
        {1,2},
        {2,3},
        {3,0},
        {1,4},
        {2,4},
        {3,4},
        {0,4}
    };
    for(i=0; i< 8; i++) {
        for(ii=0; ii< 2; ii++) {
            mPyramidEdgeIndex[i][ii]= pyedge[i][ii];
        };
    };
    uiint quedge[4][2]= {
        {0,1},
        {1,2},
        {2,3},
        {3,0}
    };
    for(i=0; i< 4; i++) {
        for(ii=0; ii< 2; ii++) {
            mQuadEdgeIndex[i][ii]= quedge[i][ii];
        };
    };
    uiint triedge[3][2]= {
        {0,1},
        {1,2},
        {2,0}
    };
    for(i=0; i< 3; i++) {
        for(ii=0; ii< 2; ii++) {
            mTriangleEdgeIndex[i][ii]= triedge[i][ii];
        };
    };
    uiint bmedge[1][2]= {
        {0,1}
    };
    for(i=0; i< 1; i++) {
        for(ii=0; ii< 2; ii++) {
            mBeamEdgeIndex[i][ii]= bmedge[i][ii];
        };
    };
    uiint hexconn[8][3]= {
        {0,3,8},
        {0,1,9},
        {1,2,10},
        {2,3,11},
        {4,7,8},
        {4,5,9},
        {5,6,10},
        {7,6,11}
    };
    for(i=0; i< 8; i++) {
        for(ii=0; ii< 3; ii++) {
            mHexaConnEdge[i][ii]= hexconn[i][ii];
        }
    }
    uiint tetconn[4][3]= {
        {2,0,3},
        {0,1,4},
        {2,1,5},
        {3,4,5}
    };
    for(i=0; i< 4; i++) {
        for(ii=0; ii< 3; ii++) {
            mTetraConnEdge[i][ii]= tetconn[i][ii];
        }
    }
    uiint priconn[6][3]= {
        {1,0,3},
        {0,2,4},
        {1,2,5},
        {3,6,8},
        {6,4,7},
        {7,8,5}
    };
    for(i=0; i< 6; i++) {
        for(ii=0; ii< 3; ii++) {
            mPrismConnEdge[i][ii]= priconn[i][ii];
        }
    }
    uiint pyconn[5][4]= {
        {0,3,7,0},
        {0,1,4,0},
        {1,2,5,1},
        {2,3,6,2},
        {4,5,6,7}
    };
    for(i=0; i< 5; i++) {
        for(ii=0; ii< 4; ii++) {
            mPyramidConnEdge[i][ii]= pyconn[i][ii];
        }
    }
    uiint quconn[4][2]= {
        {0,3},
        {0,1},
        {1,2},
        {2,3}
    };
    for(i=0; i< 4; i++) {
        for(ii=0; ii< 2; ii++) {
            mQuadConnEdge[i][ii]= quconn[i][ii];
        }
    }
    uiint triconn[3][2]= {
        {0,2},
        {0,1},
        {1,2}
    };
    for(i=0; i< 3; i++) {
        for(ii=0; ii< 2; ii++) {
            mTriangleConnEdge[i][ii]= triconn[i][ii];
        }
    }
    uiint beconn[2][1]= {
        {0},
        {0}
    };
    for(i=0; i< 2; i++) {
        for(ii=0; ii< 1; ii++) {
            mBeamConnEdge[i][ii]= beconn[i][ii];
        }
    }

    mDummyEdgeIndex= 0;//ダミー値:Pointで利用
}
CEdgeTree::~CEdgeTree()
{
}
uiint& CEdgeTree::getDummyEdgeIndex()
{
    return mDummyEdgeIndex;
}
uiint& CEdgeTree::getHexaEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    CHecMPI *pMPI=CHecMPI::Instance();

    //ofstream ofs("HexaEdgeIndex_Error.txt", ios::app);

    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(3):
            mEdgeNum=3;
            break;
        case(4):
            mEdgeNum=8;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 0");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 0" << endl;
            }
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=1;
            break;
        case(5):
            mEdgeNum=9;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 1");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 1" << endl;
            }
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(1):
            mEdgeNum=1;
            break;
        case(3):
            mEdgeNum=2;
            break;
        case(6):
            mEdgeNum=10;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 2");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 2" << endl;
            }
            break;
        }
        break;
    case(3):
        switch(localNum1) {
        case(0):
            mEdgeNum=3;
            break;
        case(2):
            mEdgeNum=2;
            break;
        case(7):
            mEdgeNum=11;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 3");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 3" << endl;
            }
            break;
        }
        break;
    case(4):
        switch(localNum1) {
        case(0):
            mEdgeNum=8;
            break;
        case(5):
            mEdgeNum=4;
            break;
        case(7):
            mEdgeNum=7;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 4");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 4" << endl;
            }
            break;
        }
        break;
    case(5):
        switch(localNum1) {
        case(4):
            mEdgeNum=4;
            break;
        case(6):
            mEdgeNum=5;
            break;
        case(1):
            mEdgeNum=9;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 5");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 5" << endl;
            }
            break;
        }
        break;
    case(6):
        switch(localNum1) {
        case(2):
            mEdgeNum=10;
            break;
        case(5):
            mEdgeNum=5;
            break;
        case(7):
            mEdgeNum=6;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 6");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 6" << endl;
            }
            break;
        }
        break;
    case(7):
        switch(localNum1) {
        case(4):
            mEdgeNum=7;
            break;
        case(6):
            mEdgeNum=6;
            break;
        case(3):
            mEdgeNum=11;
            break;
        default:
            mEdgeNum=100;
            if( pMPI->getRank()==0 ) {
                //mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalid vertex number, case 7");
                //cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                //ofs << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << " case 7" << endl;
            }
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }

    //ofs.close();

    return mEdgeNum;
}
uiint& CEdgeTree::getTetraEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    CHecMPI *pMPI=CHecMPI::Instance();

    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=2;
            break;
        case(3):
            mEdgeNum=3;
            break;
        default:
            if( pMPI->getRank()==0 ) {
                mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getTetraEdgeIndex, invalid vertex number, case 0");
                cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
            }
            mEdgeNum=100;
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=1;
            break;
        case(3):
            mEdgeNum=4;
            break;
        default:
            if( pMPI->getRank()==0 ) {
                mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getTetraEdgeIndex, invalid vertex number, case 1");
                cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
            }
            mEdgeNum=100;
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(0):
            mEdgeNum=2;
            break;
        case(1):
            mEdgeNum=1;
            break;
        case(3):
            mEdgeNum=5;
            break;
        default:
            if( pMPI->getRank()==0 ) {
                mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getTetraEdgeIndex, invalid vertex number, case 2");
                cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
            }
            mEdgeNum=100;
            break;
        }
        break;
    case(3):
        switch(localNum1) {
        case(0):
            mEdgeNum=3;
            break;
        case(1):
            mEdgeNum=4;
            break;
        case(2):
            mEdgeNum=5;
            break;
        default:
            if( pMPI->getRank()==0 ) {
                mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getTetraEdgeIndex, invalid vertex number, case 3");
                cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
            }
            mEdgeNum=100;
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getPrismEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=1;
            break;
        case(3):
            mEdgeNum=3;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=2;
            break;
        case(4):
            mEdgeNum=4;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(0):
            mEdgeNum=1;
            break;
        case(1):
            mEdgeNum=2;
            break;
        case(5):
            mEdgeNum=5;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(3):
        switch(localNum1) {
        case(0):
            mEdgeNum=3;
            break;
        case(4):
            mEdgeNum=6;
            break;
        case(5):
            mEdgeNum=8;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(4):
        switch(localNum1) {
        case(3):
            mEdgeNum=6;
            break;
        case(1):
            mEdgeNum=4;
            break;
        case(5):
            mEdgeNum=7;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(5):
        switch(localNum1) {
        case(2):
            mEdgeNum=5;
            break;
        case(3):
            mEdgeNum=8;
            break;
        case(4):
            mEdgeNum=7;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getPyramidEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(3):
            mEdgeNum=3;
            break;
        case(4):
            mEdgeNum=7;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(4):
            mEdgeNum=4;
            break;
        case(2):
            mEdgeNum=1;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(1):
            mEdgeNum=1;
            break;
        case(4):
            mEdgeNum=5;
            break;
        case(3):
            mEdgeNum=2;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(3):
        switch(localNum1) {
        case(0):
            mEdgeNum=3;
            break;
        case(2):
            mEdgeNum=2;
            break;
        case(4):
            mEdgeNum=6;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(4):
        switch(localNum1) {
        case(0):
            mEdgeNum=7;
            break;
        case(1):
            mEdgeNum=4;
            break;
        case(2):
            mEdgeNum=5;
            break;
        case(3):
            mEdgeNum=6;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getQuadEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(3):
            mEdgeNum=3;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=1;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(1):
            mEdgeNum=1;
            break;
        case(3):
            mEdgeNum=2;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(3):
        switch(localNum1) {
        case(2):
            mEdgeNum=2;
            break;
        case(0):
            mEdgeNum=3;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getTriangleEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    switch(localNum0) {
    case(0):
        switch(localNum1) {
        case(1):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=2;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(1):
        switch(localNum1) {
        case(0):
            mEdgeNum=0;
            break;
        case(2):
            mEdgeNum=1;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    case(2):
        switch(localNum1) {
        case(0):
            mEdgeNum=2;
            break;
        case(1):
            mEdgeNum=1;
            break;
        default:
            mEdgeNum=100;
            break;
        }
        break;
    default:
        mEdgeNum=200;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getBeamEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    switch(localNum0) {
    case(0):
        if(localNum1==1) {
            mEdgeNum=0;
        } else {
        }
        break;
    case(1):
        if(localNum1==0) {
            mEdgeNum=0;
        } else {
        }
        break;
    default:
        mEdgeNum=100;
        break;
    }
    return mEdgeNum;
}
uiint& CEdgeTree::getDisagTypeEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    mEdgeNum=999;
    return mEdgeNum;
}
uiint* CEdgeTree::getLocalNodeNum(const uiint& elemType, const uiint& iedge)
{
    switch(elemType) {
    case(ElementType::Hexa):
    case(ElementType::Hexa2):
        return mHexaEdgeIndex[iedge];
    case(ElementType::Tetra):
    case(ElementType::Tetra2):
        return mTetraEdgeIndex[iedge];
    case(ElementType::Prism):
    case(ElementType::Prism2):
        return mPrismEdgeIndex[iedge];
    case(ElementType::Quad):
    case(ElementType::Quad2):
        return mQuadEdgeIndex[iedge];
    case(ElementType::Triangle):
    case(ElementType::Triangle2):
        return mTriangleEdgeIndex[iedge];
    case(ElementType::Beam):
    case(ElementType::Beam2):
        return mBeamEdgeIndex[iedge];
    default:
        break;
    }
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CEdgeTree::getLocalNodeNum");
    return &pLogger->getUDummyValue();
}
