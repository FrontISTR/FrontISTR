/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommPyramid.cxx
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
#include "CommPyramid.h"
using namespace pmw;
CCommPyramid::CCommPyramid()
{
    mvNodeRank.resize(5);
    mvEdgeRank.resize(8);
    mvFaceRank.resize(5);
    mvbSend.resize(5);
    mvbRecv.resize(5);
    mvbOther.resize(5);
    mvbNodeIXCheck.resize(5);
    mvbDNodeMarking.resize(5);
    uint i;
    for(i=0; i< 5; i++){
        mvbNodeIXCheck[i]=false;
        mvbDNodeMarking[i]=false;
        mvbSend[i]=false;
        mvbRecv[i]=false;
        mvbOther[i]=false;
    }
    mvvAggCommElem.resize(5);
    mvvNeibCommElemVert.resize(5);
    mvCommNodeIndex.resize(5);
}
CCommPyramid::~CCommPyramid()
{
    ;
}
bool CCommPyramid::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Pyramid) bCoin=true;
    return bCoin;
}
