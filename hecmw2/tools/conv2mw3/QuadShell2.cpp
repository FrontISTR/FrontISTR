/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   QuadShell2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "QuadShell2.h"


CQuadShell2::CQuadShell2(void)
{
	mvNode.resize(8);
}

CQuadShell2::~CQuadShell2(void)
{
}


//vector<CNode*> CQuadShell2::getFistrFaceNode(size_t nface)
//{ 
//	//MW3Ç≈ÇÕEdgeî‘çÜ
//	vector<CNode*> vNode;
//
//	switch(nface){
//	case(1):
//		vNode.push_back(mvNode[0]);
//		vNode.push_back(mvNode[4]);
//		vNode.push_back(mvNode[1]);
//		break;
//	case(2):
//		vNode.push_back(mvNode[1]);
//		vNode.push_back(mvNode[5]);
//		vNode.push_back(mvNode[2]);
//		break;
//	case(3):
//		vNode.push_back(mvNode[2]);
//		vNode.push_back(mvNode[6]);
//		vNode.push_back(mvNode[3]);
//		break;
//	case(4):
//		vNode.push_back(mvNode[3]);
//		vNode.push_back(mvNode[7]);
//		vNode.push_back(mvNode[0]);
//		break;
//	}
//
//	return vNode;
//}
//size_t CQuadShell2::getMW3FaceNum(size_t fistr_nface)
//{ 
//	//MW3Ç≈ÇÕEdgeî‘çÜ
//	size_t mw_nedge = fistr_nface-1;
//	return mw_nedge;
//}
//
//size_t CQuadShell2::getNodeID_Fistr2MW3(size_t index)
//{
//	return CElement::getNodeID(index);
//}
//
//size_t CQuadShell2::getFaceNodeNum(size_t fistr_nface)
//{
//	return CQuad2::getFaceNodeNum(fistr_nface);
//}

string CQuadShell2::getFaceType(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::QuadShell2();
}









