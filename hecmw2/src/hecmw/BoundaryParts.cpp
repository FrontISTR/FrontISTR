/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BoundaryParts.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
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
void CBoundaryParts::resizeBNode(const uiint& res_size)
{
    mvBNode.resize(res_size);
}
void CBoundaryParts::setBNode(const uiint& ivert, CBoundaryNode* pBNode)
{
    mvBNode[ivert]= pBNode;
    mmBNodeID2Index[pBNode->getID()] = ivert;
}
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
