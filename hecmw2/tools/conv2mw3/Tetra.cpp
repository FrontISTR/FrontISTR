/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Tetra.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Tetra.h"


CTetra::CTetra(void)
{
	mvNode.resize(4);
}

CTetra::~CTetra(void)
{
}

vector<CNode*> CTetra::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		break;
	case(2):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[3]);
		break;
	case(3):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		break;
	case(4):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		break;
	}

	return vNode;
}
size_t CTetra::getMW3FaceNum(size_t fistr_nface)
{ 
	size_t mw_nface = fistr_nface-1;
	return mw_nface;
}

size_t CTetra::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}

size_t CTetra::getFaceNodeNum(size_t fistr_nface)
{
	if(fistr_nface > 4){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return 0;
	}

	return 3;
}
string CTetra::getFaceType(size_t fistr_nface)
{
	if(fistr_nface > 4){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::Triangle();
}

vector<CNode*> CTetra::getFistrEdgeNode(size_t nedge)
{
	vector<CNode*> vNode;
	//--
	//FrontISTR‚Í"1"‚©‚çŽn‚Ü‚é
	//--
	switch(nedge){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		return vNode;

	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[3]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		return vNode;

	default:
		return vNode;
	}
}


string CTetra::getEdgeType(size_t nedge)
{
	if(nedge <= 6){
		return FistrElemTypeS::Beam();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CTetra::getMW3EdgeNum(size_t fistr_nedge)
{
	if(fistr_nedge <= 6 && fistr_nedge > 0){
		return fistr_nedge-1;
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);
		return 0;
	}
}




