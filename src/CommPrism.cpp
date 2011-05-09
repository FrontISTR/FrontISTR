/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommPrism.cxx
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
#include "CommPrism.h"
using namespace pmw;
CCommPrism::CCommPrism()
{
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
    mvvAggCommElem.resize(6);
    mvvNeibCommElemVert.resize(6);
    mvCommNodeIndex.resize(6);
}
CCommPrism::~CCommPrism()
{
    ;
}
bool CCommPrism::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Prism) bCoin=true;
    return bCoin;
}
