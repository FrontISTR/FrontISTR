//
//  BoundaryParts.cpp
//
//
//          2010.06.09
//          k.Takeda
#include "BoundaryParts.h"
using namespace pmw;

CBoundaryParts::CBoundaryParts()
{
    ;
}
CBoundaryParts::~CBoundaryParts()
{
    ;
}

// BNodeのセット
// ----
void CBoundaryParts::resizeBNode(const uiint& res_size)
{
    mvBNode.resize(res_size);
}
void CBoundaryParts::setBNode(const uiint& ivert, CBoundaryNode* pBNode)
{
    mvBNode[ivert]= pBNode;

    mmBNodeID2Index[pBNode->getID()] = ivert;
}

// 頂点に(Edge,Face,Volume)IDの集合をセット
// ----
void CBoundaryParts::setupVertexElemID()
{
    CBoundaryNode* pBNode;

    uiint numOfBNode= mvBNode.size();
    uiint ibnode;

    for(ibnode=0; ibnode < numOfBNode; ibnode++){
        pBNode= mvBNode[ibnode];

        pBNode->setAggElemID(mnID);
    };
}







