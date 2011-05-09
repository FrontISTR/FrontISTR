
#include <vector>

//
//  CommPyramid.cpp
//
//
//
//                  2009.09.09
//                  2009.09.09
//                  k.Takeda
#include "CommPyramid.h"
using namespace pmw;


// construct & destruct
//
CCommPyramid::CCommPyramid()
{
    //// prolongation-Comm要素
    //mvProgCommElem.reserve(8);

    // Node rank
    mvNodeRank.resize(5);
    mvEdgeRank.resize(8);
    mvFaceRank.resize(5);

    mvbSend.resize(5);
    mvbRecv.resize(5);
    mvbOther.resize(5);

    mvbNodeIXCheck.resize(5);
    mvbDNodeMarking.resize(5);
    uint i;
    for(i=0; i< 5; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;

        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }

    // 頂点別の要素集合
    mvvAggCommElem.resize(5);
    mvvNeibCommElemVert.resize(5);

    // CommMesh内のグローバルIndex
    mvCommNodeIndex.resize(5);
}

CCommPyramid::~CCommPyramid()
{
    ;
}

// debug method :所有しているElementの型が一致するか？
//
bool CCommPyramid::isTypeCoincidence()
{
    bool bCoin(false);

    if(mpElement->getType()==ElementType::Pyramid) bCoin=true;

    return bCoin;
}




















