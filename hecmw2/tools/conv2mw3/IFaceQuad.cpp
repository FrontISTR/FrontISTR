/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   IFaceQuad.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "IFaceQuad.h"


CIFaceQuad::CIFaceQuad(void)
{
	mvNode.resize(4);
}


CIFaceQuad::~CIFaceQuad(void)
{
}

//vector<CNode*> CIFaceQuad::getFistrFaceNode(size_t nface)
//{ 
//	//MW3Ç≈ÇÕEdgeî‘çÜ
//	vector<CNode*> vNode;
//
//	switch(nface){
//	case(1):
//		vNode.push_back(mvNode[0]);
//		vNode.push_back(mvNode[1]);
//		break;
//	case(2):
//		vNode.push_back(mvNode[1]);
//		vNode.push_back(mvNode[2]);
//		break;
//	case(3):
//		vNode.push_back(mvNode[2]);
//		vNode.push_back(mvNode[3]);
//		break;
//	case(4):
//		vNode.push_back(mvNode[3]);
//		vNode.push_back(mvNode[0]);
//		break;
//	}
//
//	return vNode;
//}
//size_t CIFaceQuad::getMW3FaceNum(size_t fistr_nface)
//{ 
//	//MW3Ç≈ÇÕEdgeî‘çÜ
//	size_t mw_nedge = fistr_nface-1;
//	return mw_nedge;
//}
//
//size_t CIFaceQuad::getNodeID_Fistr2MW3(size_t index)
//{
//	return CElement::getNodeID(index);
//}
//
//size_t CIFaceQuad::getFaceNodeNum(size_t fistr_nface)
//{
//	return CQuad::getFaceNodeNum(fistr_nface);
//}

string CIFaceQuad::getFaceType(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::IFaceQuad();
}





