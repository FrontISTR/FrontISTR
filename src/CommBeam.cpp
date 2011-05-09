/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommBeam.cxx
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
#include "CommBeam.h"
using namespace pmw;
CCommBeam::CCommBeam()
{
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
