//
//  FaceTree.cpp
//
// 局所ノード番号から、面番号を提供
//
//                      2009.06.25
//                      2009.06.25
//                      k.Takeda
#include "FaceTree.h"
#include "Logger.h"
using namespace pmw;

CFaceTree::CFaceTree()
{
    // Hexa
    // --
    uint hface[6][4]={
    {0,1,2,3},
    {4,7,6,5},
    {1,5,6,2},
    {7,4,0,3},
    {0,4,5,1},
    {2,6,7,3}
    };

    uint i,ii;
    for(i=0; i< 6; i++){
    for(ii=0; ii< 4; ii++){
        mHexaFaceIndex[i][ii]= hface[i][ii];
    };
    };

    // Tetra
    // --
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

    // Prism
    // --
    uint pface[5][4]={
    {0,1,2,0},// Tri:三角形(最初のノードを最後に追加して配列要素数を合わせてある)
    {3,5,4,3},// Tri:三角形
    {0,3,4,1},// Quad:四辺形
    {1,2,5,4},// Quad:四辺形
    {0,2,5,3} // Quad:四辺形
    };

    for(i=0; i< 5; i++){
    for(ii=0; ii< 4; ii++){
       mPrismFaceIndex[i][ii]= pface[i][ii];
    };
    };

    // Pyramid
    // --
    uint pyface[5][4]={
    {0,1,2,3},// Quad:四辺形
    {1,2,4,1},// Tri:三角形(最初のノードを最後に追加して配列要素数を合わせてある)
    {2,4,3,2},// Tri:三角形
    {3,4,0,3},// Tri:三角形
    {0,4,1,0} // Tri:三角形
    };

    for(i=0; i< 5; i++){
    for(ii=0; ii< 4; ii++){
       mPyramidFaceIndex[i][ii]= pyface[i][ii];
    };
    };


    // Quad
    // --
    uint quad[4]={0,1,2,3};
    for(i=0; i< 4; i++){
        mQuadFaceIndex[i]= quad[i];
    };

    // Triangle
    // --
    uint tri[3]={0,1,2};
    for(i=0; i< 3; i++){
        mTriangleFaceIndex[i]= tri[i];
    };

}
CFaceTree::~CFaceTree()
{
//    //debug
//    cout << "~CFaceTree" << endl;
}

// 局所ノード番号(vLocalNodeIndex)から,Face(Surface)の局所インデックス番号を出力
//
// Hexa
// --
uint& CFaceTree::getHexaFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;

    // ------------------------------------------------------
    // ! 局所ノード番号がループ順序で出てくることを期待していない処理
    // ------------------------------------------------------

    // 面の数==6, 面の頂点数==4
    uint ia,ib;
    for(faceIndex=0; faceIndex< 6; faceIndex++){
        hitCount=0;//初期化
        for(ia=0; ia< 4; ia++){
        for(ib=0; ib< 4; ib++){
            if(mHexaFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]){ ++hitCount; break;}
        };
        };
        //if(hitCount==4) return faceIndex;//頂点の数だけ一致
        if(hitCount==4){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
// Hexa
// --
// 節点のツリーによる分岐 -> Faceを特定(四辺形なので,起点と対角点で面を特定)
//
uint& CFaceTree::getHexaFaceIndex2(const vuint& vLocalNodeIndex)
{
    switch(vLocalNodeIndex[0]){
    case(0)://起点 0
        switch(vLocalNodeIndex[2]){//対角点
            case(2): mFaceIndex=0; break;
            case(5): mFaceIndex=4; break;
            case(7): mFaceIndex=3; break;           
        }
        break;
    case(1):// 起点 1
        switch(vLocalNodeIndex[2]){//対角点
            case(3): mFaceIndex=0; break;
            case(6): mFaceIndex=2; break;
            case(4): mFaceIndex=4; break;
        }
        break;
    case(2):// 起点 2
        switch(vLocalNodeIndex[2]){//対角点
            case(5): mFaceIndex=2; break;
            case(0): mFaceIndex=0; break;
            case(7): mFaceIndex=5; break;
        }
        break;
    case(3):// 起点 3
        switch(vLocalNodeIndex[2]){//対角点
            case(1): mFaceIndex=0; break;
            case(6): mFaceIndex=5; break;
            case(4): mFaceIndex=3; break;
        }
        break;
    case(4):// 起点 4
        switch(vLocalNodeIndex[2]){//対角点
            case(1): mFaceIndex=4; break;
            case(3): mFaceIndex=3; break;
            case(6): mFaceIndex=1; break;
        }
        break;
    case(5):// 起点 5
        switch(vLocalNodeIndex[2]){//対角点
            case(2): mFaceIndex=2; break;
            case(0): mFaceIndex=4; break;
            case(7): mFaceIndex=1; break;
        }
        break;
    case(6):// 起点 6
        switch(vLocalNodeIndex[2]){//対角点
            case(1): mFaceIndex=2; break;
            case(3): mFaceIndex=5; break;
            case(4): mFaceIndex=1; break;
        }
        break;
    case(7):// 起点 7
        switch(vLocalNodeIndex[2]){//対角点
            case(2): mFaceIndex=5; break;
            case(5): mFaceIndex=1; break;
            case(0): mFaceIndex=3; break;
        }
        break;
    }
    return mFaceIndex;
}




// Tetra
// --
uint& CFaceTree::getTetraFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;

    // 面の数==4, 面の頂点数==3
    uint ia,ib;
    for(faceIndex=0; faceIndex< 4; faceIndex++){
        hitCount=0;//初期化
        for(ia=0; ia< 3; ia++){
        for(ib=0; ib< 3; ib++){
            if(mTetraFaceIndex[faceIndex][ia]==vLocalNodeIndex[ib]) ++hitCount;
        };
        };
        //if(hitCount==3) return faceIndex;//頂点の数だけ一致した面番号
        if(hitCount==3){
            mFaceIndex=faceIndex;
            return mFaceIndex;
        }
    };
}
// Tetra
// --
// 節点のツリーによる分岐
//
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



// Prism
// --
uint& CFaceTree::getPrismFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;

    // 面の数==5, 面の頂点数==4
    //   0番め,1番めのFaceは三角形であるが、頂点数を揃えるため4頂点にしてある.
    //   三角形の面は、最初の局所ノード番号を最後に再度使ってある.
    // => 3つの頂点のみチェックする.
    uint ia,ib;
    for(faceIndex=0; faceIndex< 5; faceIndex++){
        hitCount=0;//初期化
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
// Prism
// --
// 節点のツリーによる分岐
//
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


// Pyramid
// --
uint& CFaceTree::getPyramidFaceIndex(const vuint& vLocalNodeIndex)
{
    uint faceIndex;
    uint hitCount;

    // 面の数==5, 面の頂点数==4
    //   最初以外の、4つののFaceは三角形であるが、頂点数を揃えるため4頂点にしてある.
    //   三角形の面は、最初の局所ノード番号を最後に再度使ってある. ＝＞ 面は,結局3つの頂点のみチェックするので関係ない.
    uint ia,ib;
    for(faceIndex=0; faceIndex< 5; faceIndex++){
        hitCount=0;//初期化
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

// Pyramid
// --
// 節点のツリー構造による面探索
//
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



// Quad_Shell
// --
// ノードが一致すれば,面番号=”0”を返す.
//
uint& CFaceTree::getQuadFaceIndex(const vuint& vLocalNodeIndex)
{
    // 面の頂点数==4であるが,
    //  面判定は3個の頂点で判定するので引数vuintのサイズ数を基準とする.
    uint ia,ib;
    uint hitCount=0;//初期化
    uint numOfLocalNode = vLocalNodeIndex.size();

    for(ia=0; ia< 4; ia++){
    for(ib=0; ib< numOfLocalNode; ib++){
        if(mQuadFaceIndex[ia]==vLocalNodeIndex[ib]) ++hitCount;
    };
    };
    
    // 正否-判定
    if(hitCount >= numOfLocalNode){
        mFaceIndex=0;//正しい面番号
    }else{
        mFaceIndex=1;//間違い
    }

    return mFaceIndex;
}


// Triangle_Shell
// --
// ノードが一致すれば,面番号=”0”を返す.
//
uint& CFaceTree::getTriangleFaceIndex(const vuint& vLocalNodeIndex)
{
    // 面の頂点数==3
    uint ia,ib;
    uint hitCount=0;//初期化
    uint numOfLocalNode = vLocalNodeIndex.size();
    
    for(ia=0; ia< 3; ia++){
    for(ib=0; ib< numOfLocalNode; ib++){
        if(mTriangleFaceIndex[ia]==vLocalNodeIndex[ib]) ++hitCount;
    };
    };


    // 正否-判定
    if(hitCount==numOfLocalNode){
        mFaceIndex=0;//正しい面番号
    }else{
        mFaceIndex=1;//間違い
    }

    return mFaceIndex;
}


// Faceを表す局所ノード番号
// --
// Hexa
// --
uint* CFaceTree::getLocalNodeHexaFace(const uint& faceIndex)
{
    if(faceIndex >= 6){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodeHexaFace",faceIndex);
    }

    return mHexaFaceIndex[faceIndex];
}
// Tetra
// --
uint* CFaceTree::getLocalNodeTetraFace(const uint& faceIndex)
{
    if(faceIndex >= 4){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodeTetraFace",faceIndex);
    }

    return mTetraFaceIndex[faceIndex];
}
// Prism
// --
uint* CFaceTree::getLocalNodePrismFace(const uint& faceIndex)
{
    if(faceIndex >= 5){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodePrismFace",faceIndex);
    }

    return mPrismFaceIndex[faceIndex];
}
// Pyramid
// --
uint* CFaceTree::getLocalNodePyramidFace(const uint& faceIndex)
{
    if(faceIndex >= 5){
        Utility::CLogger *pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error,"faceIndex Over @CFaceTree::getLocalNodePyramidFace",faceIndex);
    }

    return mPyramidFaceIndex[faceIndex];
}
// Quad
// --
uint* CFaceTree::getLocalNodeQuadFace()
{
    return mQuadFaceIndex;
}
// Triangle
// --
uint* CFaceTree::getLocalNodeTriangleFace()
{
    return mTriangleFaceIndex;
}



