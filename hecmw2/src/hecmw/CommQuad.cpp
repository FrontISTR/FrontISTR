/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommQuad.cpp
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
#include "CommQuad.h"
using namespace pmw;
CCommQuad::CCommQuad()
{
    mvNodeRank.resize(4);
    mvEdgeRank.resize(4);
    mvFaceRank.resize(1);
    mvbSend = new bool[4];
    mvbRecv = new bool[4];
    mvbOther = new bool[4];
    mvbNodeIXCheck = new bool[4];
    mvbDNodeMarking = new bool[4];
    uiint i;
    for(i=0; i< 4; i++) {
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
CCommQuad::~CCommQuad()
{
    ;
}
bool CCommQuad::isTypeCoincidence()
{
    bool bCoin(false);
    if(mpElement->getType()==ElementType::Quad) bCoin=true;
    return bCoin;
}
