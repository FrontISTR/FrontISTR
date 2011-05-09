/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommTetra.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <vector>
#include "CommTetra.h"
using namespace pmw;
CCommTetra::CCommTetra()
{
    mvNodeRank.resize(4);
    mvEdgeRank.resize(6);
    mvFaceRank.resize(4);
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
    mvvAggCommElem.resize(4);
    mvvNeibCommElemVert.resize(4);
    mvCommNodeIndex.resize(4);
}
CCommTetra::~CCommTetra()
{
}
bool CCommTetra::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Tetra) bCoin=true;
    return bCoin;
}
