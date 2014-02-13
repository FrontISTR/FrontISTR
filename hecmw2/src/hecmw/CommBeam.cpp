/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommBeam.cpp
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
#include "CommBeam.h"
using namespace pmw;
CCommBeam::CCommBeam()
{
    mvNodeRank.resize(2);
    mvEdgeRank.resize(1);
    mvbSend = new bool[2];
    mvbRecv = new bool[2];
    mvbOther = new bool[2];
    mvbNodeIXCheck = new bool[2];
    mvbDNodeMarking = new bool[2];
    uiint i;
    for(i=0; i< 2; i++) {
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;
        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }
    mvvAggCommElem.resize(2);
    mvvNeibCommElemVert.resize(2);
    mvCommNodeIndex.resize(2);
}
CCommBeam::~CCommBeam()
{
    ;
}
bool CCommBeam::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Beam) bCoin=true;
    return bCoin;
}
