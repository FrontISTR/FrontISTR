
#include <vector>

//
//  CommPrism.cpp
//
//
//
//                  2009.09.08
//                  2009.09.08
//                  k.Takeda
#include "CommPrism.h"
using namespace pmw;

// construct & destruct
//
CCommPrism::CCommPrism()
{
    //// prolongation-Comm要素
    //mvProgCommElem.reserve(6);

    // Node rank
    mvNodeRank.resize(6);
    mvEdgeRank.resize(9);
    mvFaceRank.resize(5);

    mvbSend.resize(6);
    mvbRecv.resize(6);
    mvbOther.resize(6);

    mvbNodeIXCheck.resize(6);
    mvbDNodeMarking.resize(6);
    uint i;
    for(i=0; i< 6; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;

        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }

    // 頂点別の要素集合
    mvvAggCommElem.resize(6);
    mvvNeibCommElemVert.resize(6);

    // CommMesh内のグローバルIndex
    mvCommNodeIndex.resize(6);
}

CCommPrism::~CCommPrism()
{
    ;
}


// debug method :所有しているElementの型が一致するか？
//
bool CCommPrism::isTypeCoincidence()
{
    bool bCoin(false);

    if(mpElement->getType()==ElementType::Prism) bCoin=true;

    return bCoin;
}












