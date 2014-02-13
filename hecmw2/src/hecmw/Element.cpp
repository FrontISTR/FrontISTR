/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Element.cpp
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
#include "Element.h"
using namespace pmw;
CElement::CElement(void)
{
    mbComm =false;
    mbDComm=false;
    mbRComm=false;
    mCommID= -1;
    mbMaster=false;
    mbSlave= false;
    mbCommMesh2=false;
}
CElement::~CElement(void)
{
}
void CElement::setNode(CNode* pNode, const uiint& local_id)
{
    mvNode[local_id] = pNode;
    mmIDLocal[pNode->getID()]= local_id;
}

void CElement::ReInit_IDLocal()
{
    uiint i, nNumOfNode = mvNode.size();
    for(i=0; i < nNumOfNode; i++) {
        uiint id = mvNode[i]->getID();
        mmIDLocal[id] = i;
    };
}
void CElement::reserveEdgeElement(const uiint& edgeIndex, const uiint& numOfElem)
{
    string s_msg("reserveEdgeElement");
    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)) {
        mvvEdgeElement[edgeIndex].reserve(numOfElem);
    } else {
        return;
    }
}
void CElement::setEdgeElement(const uiint& edgeIndex, CElement* pElem)
{
    string s_msg("setEdgeElement");
    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)) {
        mvvEdgeElement[edgeIndex].push_back(pElem);
    } else {
        return;
    }
}
void CElement::setEdgeAggElement(const uiint& edgeIndex, vector<CElement*> vElement)
{
    mvvEdgeElement[edgeIndex] = vElement;
}
void CElement::setEdgeInterNode(CNode* pNode, const uiint& edgeIndex)
{
    mvEdgeInterNode[edgeIndex]=pNode;
}
vuint CElement::getLocalNodeNum(CNode* pNode0, CNode* pNode1)
{
    uiint id0, id1;
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    vuint pairIndex;
    pairIndex.resize(2);
    pairIndex[0]= mmIDLocal[id0];
    pairIndex[1]= mmIDLocal[id1];
    return pairIndex;
}
bool CElement::isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uiint iface= getFaceIndex(pNode0,pNode1,pNode2);
    return mvb_face[iface];
}
void CElement::setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uiint iface = getFaceIndex(pNode0,pNode1,pNode2);
    mvb_face[iface]=true;
}
CElement* CElement::getProgElem_NodeID(const uiint& nodeID)
{
    uiint ivert= mmIDLocal[nodeID];
    return mvProgElement[ivert];
}


void CElement::markingCommFace(const uiint& iface)
{
    mvCommLFace.push_back(iface);
    mvbCommLFace[iface]=true;
}
void CElement::markingCommEdge(const uiint& iedge)
{
    mvCommLEdge.push_back(iedge);
    mvbCommLEdge[iedge]=true;
}
void CElement::markingCommVert(const uiint& ivert)
{
    mvCommLVert.push_back(ivert);
    mvbCommLVert[ivert]=true;
}





