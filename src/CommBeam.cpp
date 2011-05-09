
#include <vector>

//
//  CommBeam.cpp
//
//
//
//                  2009.09.10
//                  2009.09.10
//                  k.Takeda
#include "CommBeam.h"
using namespace pmw;

// construct & destruct
//
CCommBeam::CCommBeam()
{
    //// prolongation-Comm要素
    //mvProgCommElem.reserve(2);

    // Node rank
    mvNodeRank.resize(2);
    mvEdgeRank.resize(1);

    mvbSend.resize(2);
    mvbRecv.resize(2);
    mvbOther.resize(2);

    mvbNodeIXCheck.resize(2);
    mvbDNodeMarking.resize(2);
    uint i;
    for(i=0; i< 2; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;

        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }

    // 頂点別の要素集合
    mvvAggCommElem.resize(2);
    mvvNeibCommElemVert.resize(2);

    // CommMesh内のグローバルIndex
    mvCommNodeIndex.resize(2);
}

CCommBeam::~CCommBeam()
{
    ;
}

// debug method :所有しているElementの型が一致するか？
//
bool CCommBeam::isTypeCoincidence()
{
    bool bCoin(false);

    if(mpElement->getType()==ElementType::Beam) bCoin=true;

    return bCoin;
}



























