/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommTriangle.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
    mvbSend = new bool[3];
    mvbRecv = new bool[3];
    mvbOther = new bool[3];
    mvbNodeIXCheck = new bool[3];
    mvbDNodeMarking = new bool[3];
    uiint i;
    for(i=0; i< 3; i++) {
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
