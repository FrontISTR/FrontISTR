/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommTriangle.cxx
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
#include "CommTriangle.h"
using namespace pmw;
CCommTriangle::CCommTriangle()
{
    mvNodeRank.resize(3);
    mvEdgeRank.resize(3);
    mvFaceRank.resize(1);
    mvbSend.resize(3);
    mvbRecv.resize(3);
    mvbOther.resize(3);
    mvbNodeIXCheck.resize(3);
    mvbDNodeMarking.resize(3);
    uint i;
    for(i=0; i< 3; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;
        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }
    mvvAggCommElem.resize(3);
    mvvNeibCommElemVert.resize(3);
    mvCommNodeIndex.resize(3);
}
CCommTriangle::~CCommTriangle()
{
    ;
}
bool CCommTriangle::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Triangle) bCoin=true;
    return bCoin;
}
