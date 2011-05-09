
#include <vector>

//
//  CommQuad.cpp
//
//
//
//                  2009.09.09
//                  2009.09.09
//                  k.Takeda
#include "CommQuad.h"
using namespace pmw;

// construct & destruct
//
CCommQuad::CCommQuad()
{
    //// prolongation-Comm要素
    //mvProgCommElem.reserve(4);

    // Node rank
    mvNodeRank.resize(4);
    mvEdgeRank.resize(4);
    mvFaceRank.resize(1);

    mvbSend.resize(4);
    mvbRecv.resize(4);
    mvbOther.resize(4);

    mvbNodeIXCheck.resize(4);
    mvbDNodeMarking.resize(4);
    uint i;
    for(i=0; i< 4; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;

        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }

    // 頂点別の要素集合
    mvvAggCommElem.resize(4);
    mvvNeibCommElemVert.resize(4);

    // CommMesh内のグローバルIndex
    mvCommNodeIndex.resize(4);
}

CCommQuad::~CCommQuad()
{
    ;
}


// debug method :所有しているElementの型が一致するか？
//
bool CCommQuad::isTypeCoincidence()
{
    bool bCoin(false);

    if(mpElement->getType()==ElementType::Quad) bCoin=true;

    return bCoin;
}















