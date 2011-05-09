/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Element.cxx
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
#include "Element.h"
using namespace pmw;
CElement::CElement(void)
{
    mvPairNodeLocalNum.resize(2);
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
void CElement::setNode(CNode* pNode, const uint& local_id)
{
    mvNode[local_id] = pNode;
    mmIDLocal[pNode->getID()]= local_id;
}
void CElement::reserveEdgeElement(const uint& edgeIndex, const uint& numOfElem)
{
    string s_msg("reserveEdgeElement");
    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)){
        mvvEdgeElement[edgeIndex].reserve(numOfElem);
    }else{
        return;
    }
}
void CElement::setEdgeElement(const uint& edgeIndex, CElement* pElem)
{
    string s_msg("setEdgeElement");
    if(IndexCheck(ElementPropType::Edge, edgeIndex, s_msg)){
        mvvEdgeElement[edgeIndex].push_back(pElem);
    }else{
        return;
    }
}
void CElement::setEdgeAggElement(const uint& edgeIndex, vector<CElement*> vElement)
{
    mvvEdgeElement[edgeIndex] = vElement;
}
void CElement::setEdgeInterNode(CNode* pNode, const uint& edgeIndex)
{
    mvEdgeInterNode[edgeIndex]=pNode;
}
vuint& CElement::getLocalNodeNum(CNode* pNode0, CNode* pNode1)
{
    uint id0, id1;
    id0 = pNode0->getID();
    id1 = pNode1->getID();
    mvPairNodeLocalNum[0]= mmIDLocal[id0];
    mvPairNodeLocalNum[1]= mmIDLocal[id1];
    return mvPairNodeLocalNum;
}
bool CElement::isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uint iface= getFaceIndex(pNode0,pNode1,pNode2);
    return mvb_face[iface];
}
void CElement::setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2)
{
    uint iface = getFaceIndex(pNode0,pNode1,pNode2);
    mvb_face[iface]=true;
}
CElement* CElement::getProgElem_NodeID(const uint& nodeID)
{
    uint ivert= mmIDLocal[nodeID];
    return mvProgElement[ivert];
}
