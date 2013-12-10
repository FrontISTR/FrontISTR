/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Triangle2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Triangle2.h"


CTriangle2::CTriangle2(void)
{
	mvNode.resize(6);
}

CTriangle2::~CTriangle2(void)
{
}

vector<CNode*> CTriangle2::getFistrFaceNode(size_t nface)
{ 
	//MW3Ç≈ÇÕEdgeî‘çÜ
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[1]);
		break;
	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[2]);
		break;
	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[0]);
		break;
	}

	return vNode;
}
size_t CTriangle2::getMW3FaceNum(size_t fistr_nface)
{ 
	//MW3Ç≈ÇÕEdgeî‘çÜ
	size_t mw_nedge = fistr_nface-1;
	return mw_nedge;
}

size_t CTriangle2::getNodeID_Fistr2MW3(size_t index)
{
	CMessage *pMsg=CMessage::Instance();

	size_t up_index;

	switch(index){
	case(0):
		up_index= 0; break;
	case(1):
		up_index= 1; break;
	case(2):
		up_index= 2; break;
	case(3):
		up_index= 4; break;
	case(4):
		up_index= 5; break;
	case(5):
		up_index= 3; break;
	default:
		string sIX = lexical_cast<string>(index);
		pMsg->error(sIX+" was out of range, Triangle2::getNodeID_Fistr2MW3");
		pMsg->info("index set to 5");

		up_index= 5;
		break;
	}

	return mvNode[up_index]->getID();
}

size_t CTriangle2::getFaceNodeNum(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceNodeNum, arg was out of range, elem id:"+sID);
		return 0;
	}

	return 6;
}
string CTriangle2::getFaceType(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::Triangle2();
}

vector<CNode*> CTriangle2::getFistrEdgeNode(size_t nedge)
{
	vector<CNode*> vNode;
	//--
	//FrontISTRÇÕ"1"Ç©ÇÁénÇ‹ÇÈ
	//--
	switch(nedge){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[4]);
		return vNode;

	default:
		return vNode;
	}
}

string CTriangle2::getEdgeType(size_t nedge)
{
	if(nedge <= 3){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CTriangle2::getMW3EdgeNum(size_t fistr_nedge)
{
	if(fistr_nedge <= 3 && fistr_nedge > 0){
		return fistr_nedge-1;
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);
		return 0;
	}
}



