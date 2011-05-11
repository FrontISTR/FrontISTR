//
//  EdgeTree.cpp
//  各種の要素形状の、辺(Edge)のノード接続
//
//
//                      2009.09.01
//                      2009.06.23
//                      k.Takeda
#include "EdgeTree.h"
#include "ElementType.h"
#include "Logger.h"
using namespace pmw;

#include <iostream>
//
//
CEdgeTree::CEdgeTree()
{
    mpLogger = Utility::CLogger::Instance();

    uiint i,ii;
    // Edgeを構成する両端の局所ノード番号
    //--
    //Hexa
    uiint hexedge[12][2]={
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
    for(i=0; i< 12; i++){
        for(ii=0; ii< 2; ii++){
            mHexaEdgeIndex[i][ii]= hexedge[i][ii];
        };
    };

    //Tetra
    uiint tetedge[6][2]={
        {0,1},
        {1,2},
        {2,0},
        {0,3},
        {1,3},
        {2,3}
    };
    for(i=0; i< 6; i++){
        for(ii=0; ii< 2; ii++){
            mTetraEdgeIndex[i][ii]= tetedge[i][ii];
        };
    };

    //Prism
    uiint priedge[9][2]={
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
    for(i=0; i< 9; i++){
        for(ii=0; ii< 2; ii++){
            mPrismEdgeIndex[i][ii]= priedge[i][ii];
        };
    };

    //Pyramid
    uiint pyedge[8][2]={
        {0,1},
        {1,2},
        {2,3},
        {3,0},
        {1,4},
        {2,4},
        {3,4},
        {0,4}
    };
    for(i=0; i< 8; i++){
        for(ii=0; ii< 2; ii++){
            mPyramidEdgeIndex[i][ii]= pyedge[i][ii];
        };
    };

    //Quad
    uiint quedge[4][2]={
        {0,1},
        {1,2},
        {2,3},
        {3,0}
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 2; ii++){
            mQuadEdgeIndex[i][ii]= quedge[i][ii];
        };
    };

    //Triangle
    uiint triedge[3][2]={
        {0,1},
        {1,2},
        {2,0}
    };
    for(i=0; i< 3; i++){
        for(ii=0; ii< 2; ii++){
            mTriangleEdgeIndex[i][ii]= triedge[i][ii];
        };
    };

    //Beam
    uiint bmedge[1][2]={
        {0,1}
    };
    for(i=0; i< 1; i++){
        for(ii=0; ii< 2; ii++){
            mBeamEdgeIndex[i][ii]= bmedge[i][ii];
        };
    };


    // ノードに接続するEdge番号
    // --
    // Hexa
    uiint hexconn[8][3]={
        {0,3,8},
        {0,1,9},
        {1,2,10},
        {2,3,11},
        {4,7,8},
        {4,5,9},
        {5,6,10},
        {7,6,11}
    };
    for(i=0; i< 8; i++){
        for(ii=0; ii< 3; ii++){
            mHexaConnEdge[i][ii]= hexconn[i][ii];
        }
    }
    
    // Tetra
    uiint tetconn[4][3]={
        {2,0,3},
        {0,1,4},
        {2,1,5},
        {3,4,5}
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 3; ii++){
            mTetraConnEdge[i][ii]= tetconn[i][ii];
        }
    }
    
    // Prism
    uiint priconn[6][3]={
        {1,0,3},
        {0,2,4},
        {1,2,5},
        {3,6,8},
        {6,4,7},
        {7,8,5}
    };
    for(i=0; i< 6; i++){
        for(ii=0; ii< 3; ii++){
            mPrismConnEdge[i][ii]= priconn[i][ii];
        }
    }
    
    // Pyramid
    uiint pyconn[5][4]={
        {0,3,7,0},
        {0,1,4,0},
        {1,2,5,1},
        {2,3,6,2},
        {4,5,6,7}
    };
    for(i=0; i< 5; i++){
        for(ii=0; ii< 4; ii++){
            mPyramidConnEdge[i][ii]= pyconn[i][ii];
        }
    }
    
    // Quad
    uiint quconn[4][2]={
        {0,3},
        {0,1},
        {1,2},
        {2,3}
    };
    for(i=0; i< 4; i++){
        for(ii=0; ii< 2; ii++){
            mQuadConnEdge[i][ii]= quconn[i][ii];
        }
    }
    
    // Triangle
    uiint triconn[3][2]={
        {0,2},
        {0,1},
        {1,2}
    };
    for(i=0; i< 3; i++){
        for(ii=0; ii< 2; ii++){
            mTriangleConnEdge[i][ii]= triconn[i][ii];
        }
    }
    
    // Beam
    uiint beconn[2][1]={
        {0},
        {0}
    };
    for(i=0; i< 2; i++){
        for(ii=0; ii< 1; ii++){
            mBeamConnEdge[i][ii]= beconn[i][ii];
        }
    }
    
}
CEdgeTree::~CEdgeTree()
{
    ////debug
    //std::cout << "~CEdgeTree" << std::endl;
}

// Hexa 辺のノード接続
//
uiint& CEdgeTree::getHexaEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,3,4 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 0");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,2,5 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 1");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(2):
            // 2番 -> 1,3,6 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 2");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(3):
            // 3番 -> 0,2,7 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 3");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(4):
            // 4番 -> 0,5,7 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 4");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(5):
            // 5番 -> 4,6,1 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 5");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(6):
            // 6番 -> 2,5,7 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 6");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        case(7):
            // 7番 -> 4,6,3 の三本の辺
            switch(localNum1){
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
                    mpLogger->Info(Utility::LoggerMode::Error, "EdgeTree::getHexaEdgeIndex, invalit vertex number, case 7");
                    cout << "EdgeTree::getHexaEdgeIndex, Num0 " << localNum0 << ",  Num1 " << localNum1 << endl;
                    break;
            }
            break;
        default:
            break;
    }
    return mEdgeNum;
}

// Tetra 辺のノード接続
//
uiint& CEdgeTree::getTetraEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,2,3 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,2,3 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(2):
            // 2番 -> 0,1,3 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(3):
            // 3番 -> 0,1,2 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;

        default:
            break;
    }
    return mEdgeNum;
}

// Prism 辺のノード接続
//
uiint& CEdgeTree::getPrismEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,2,3 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,2,4 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(2):
            // 2番 -> 0,1,5 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(3):
            // 3番 -> 0,4,5 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(4):
            // 4番 -> 3,1,5 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(5):
            // 5番 -> 2,3,4 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;

        default:
            break;
    }
    return mEdgeNum;
}

// Pyramid 辺のノード接続
//
uiint& CEdgeTree::getPyramidEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,3,4 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,4,2 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(2):
            // 2番 -> 1,4,3 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(3):
            // 3番 -> 0,2,4 の三本の辺
            switch(localNum1){
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
                    break;
            }
            break;
        case(4):
            // 4番 -> 0,1,2,3 の四本の辺
            switch(localNum1){
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
                    break;
            }
            break;

        default:
            break;
    }
    return mEdgeNum;
}


// Quad 辺のノード接続
//
uiint& CEdgeTree::getQuadEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,3 の2本の辺
            switch(localNum1){
                case(1):
                    mEdgeNum=0;
                    break;
                case(3):
                    mEdgeNum=3;
                    break;
                default:
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,2の三本の辺
            switch(localNum1){
                case(0):
                    mEdgeNum=0;
                    break;
                case(2):
                    mEdgeNum=1;
                    break;
                default:
                    break;
            }
            break;
        case(2):
            // 2番 -> 1,3 の三本の辺
            switch(localNum1){
                case(1):
                    mEdgeNum=1;
                    break;
                case(3):
                    mEdgeNum=2;
                    break;
                default:
                    break;
            }
            break;
        case(3):
            // 3番 -> 2,0の三本の辺
            switch(localNum1){
                case(2):
                    mEdgeNum=2;
                    break;
                case(0):
                    mEdgeNum=3;
                    break;
                default:
                    break;
            }
            break;

        default:
            break;
    }
    return mEdgeNum;
}

// Triangle 辺のノード接続
//
uiint& CEdgeTree::getTriangleEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1,2 の2本の辺
            switch(localNum1){
                case(1):
                    mEdgeNum=0;
                    break;
                case(2):
                    mEdgeNum=2;
                    break;
                default:
                    break;
            }
            break;
        case(1):
            // 1番 -> 0,2の三本の辺
            switch(localNum1){
                case(0):
                    mEdgeNum=0;
                    break;
                case(2):
                    mEdgeNum=1;
                    break;
                default:
                    //Error
                    break;
            }
            break;
        case(2):
            // 2番 -> 0,1の三本の辺
            switch(localNum1){
                case(0):
                    mEdgeNum=2;
                    break;
                case(1):
                    mEdgeNum=1;
                    break;
                default:
                    //Error
                    break;
            }
            break;
        default:
            break;
    }
    return mEdgeNum;
}

// Beam 辺のノード接続
//
uiint& CEdgeTree::getBeamEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    // Edge(辺)のツリー
    switch(localNum0){
        case(0):
            // 局所番号0番 -> 一方の端点の局所番号 1 の1本の辺
            if(localNum1==1){
                mEdgeNum=0;
            }else{
                //Error
            }
            break;
        case(1):
            // 局所番号1番 -> 一方の端点の局所番号 0 の1本の辺
            if(localNum1==0){
                mEdgeNum=0;
            }else{
                //Error
            }
            break;
        default:
            // Error
            break;
    }

    return mEdgeNum;
}

// Disagreement(不一致)Typeの場合の辺番号:存在しない辺番号
//
uiint& CEdgeTree::getDisagTypeEdgeIndex(const uiint& localNum0, const uiint& localNum1)
{
    mEdgeNum=999;
    return mEdgeNum;
}



// 辺(Edge)番号から,端点の局所ノード番号を提供: 統合版
// ----
uiint* CEdgeTree::getLocalNodeNum(const uiint& elemType, const uiint& iedge)
{
    switch(elemType){
        case(ElementType::Hexa):case(ElementType::Hexa2):

            return mHexaEdgeIndex[iedge];

        case(ElementType::Tetra):case(ElementType::Tetra2):

            return mTetraEdgeIndex[iedge];

        case(ElementType::Prism):case(ElementType::Prism2):

            return mPrismEdgeIndex[iedge];

        case(ElementType::Quad):case(ElementType::Quad2):

            return mQuadEdgeIndex[iedge];

        case(ElementType::Triangle):case(ElementType::Triangle2):

            return mTriangleEdgeIndex[iedge];

        case(ElementType::Beam):case(ElementType::Beam2):

            return mBeamEdgeIndex[iedge];

        default:
            break;
    }

    //default処理
    //----
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Error, "invalid ElementType, CEdgeTree::getLocalNodeNum");

    return &pLogger->getUDummyValue();
}

















