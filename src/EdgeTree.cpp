//
//  EdgeTree.cpp
//  各種の要素形状の、辺(Edge)のノード接続
//
//
//                      2009.06.23
//                      2009.06.23
//                      k.Takeda
#include "EdgeTree.h"
using namespace pmw;

#include <iostream>
//
//
CEdgeTree::CEdgeTree()
{
    ;
}
CEdgeTree::~CEdgeTree()
{
    ////debug
    //std::cout << "~CEdgeTree" << std::endl;
}

// Hexa 辺のノード接続
//
uint& CEdgeTree::getHexaEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getTetraEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getPrismEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getPyramidEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getQuadEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getTriangleEdgeIndex(const uint& localNum0, const uint& localNum1)
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
uint& CEdgeTree::getBeamEdgeIndex(const uint& localNum0, const uint& localNum1)
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










