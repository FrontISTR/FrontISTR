/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Quad2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Quad2.h"


CQuad2::CQuad2(void)
{
	mvNode.resize(8);
}

CQuad2::~CQuad2(void)
{
}

vector<CNode*> CQuad2::getFistrFaceNode(size_t nface)
{ 
	//MW3Ç≈ÇÕEdgeî‘çÜ
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[1]);
		break;
	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[2]);
		break;
	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[3]);
		break;
	case(4):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[0]);
		break;
	}

	return vNode;
}
size_t CQuad2::getMW3FaceNum(size_t fistr_nface)
{ 
	//MW3Ç≈ÇÕEdgeî‘çÜ
	size_t mw_nedge = fistr_nface-1;
	return mw_nedge;
}

size_t CQuad2::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}

size_t CQuad2::getFaceNodeNum(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceNodeNum, arg was out of range, elem id:"+sID);
		return 0;
	}

	return 8;
}
string CQuad2::getFaceType(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::Quad2();
}

vector<CNode*> CQuad2::getFistrEdgeNode(size_t nedge)
{
	vector<CNode*> vNode;
	//--
	//FrontISTRÇÕ"1"Ç©ÇÁénÇ‹ÇÈ
	//--
	switch(nedge){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[6]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[7]);
		return vNode;

	default:
		return vNode;
	}
}

string CQuad2::getEdgeType(size_t nedge)
{
	if(nedge <= 4){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CQuad2::getMW3EdgeNum(size_t fistr_nedge)
{
	if(fistr_nedge <= 4 && fistr_nedge > 0){
		return fistr_nedge-1;
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);
		return 0;
	}
}


