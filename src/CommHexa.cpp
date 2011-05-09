/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommHexa.cxx
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
#include "CommHexa.h"
using namespace pmw;
CCommHexa::CCommHexa()
{
    mvNodeRank.resize(8);
    mvEdgeRank.resize(12);
    mvFaceRank.resize(6);
    mvbSend.resize(8);
    mvbRecv.resize(8);
    mvbOther.resize(8);
    mvbNodeIXCheck.resize(8);
    mvbDNodeMarking.resize(8);
    uint i;
    for(i=0; i< 8; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;
        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }
    mvvAggCommElem.resize(8);
    mvvNeibCommElemVert.resize(8);
    mvCommNodeIndex.resize(8);
}
CCommHexa::~CCommHexa()
{
    ;
}
bool CCommHexa::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Hexa) bCoin=true;
    return bCoin;
}
