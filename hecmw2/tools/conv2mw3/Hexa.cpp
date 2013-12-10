/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Hexa.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Hexa.h"


CHexa::CHexa(void)
{
	mvNode.resize(8);
}

CHexa::~CHexa(void)
{
}

vector<CNode*> CHexa::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		break;
	case(2):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[7]);
		break;
	case(3):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[4]);
		break;
	case(4):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[5]);
		break;
	case(5):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[6]);
		break;
	case(6):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[7]);
		break;
	}

	return vNode;
}
size_t CHexa::getMW3FaceNum(size_t fistr_nface)
{ 
	size_t mw_nface;

	switch(fistr_nface){
	case(1):
		mw_nface= 0;
		break;
	case(2):
		mw_nface= 1;
		break;
	case(3):
		mw_nface= 4;
		break;
	case(4):
		mw_nface= 2;
		break;
	case(5):
		mw_nface= 5;
		break;
	case(6):
		mw_nface= 3;
		break;
	}

	return mw_nface;
}

size_t CHexa::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}

size_t CHexa::getFaceNodeNum(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return 4;
	case(2):
		return 4;
	case(3):
		return 4;
	case(4):
		return 4;
	case(5):
		return 4;
	case(6):
		return 4;
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return 0;
	}
}
string CHexa::getFaceType(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return FistrElemTypeS::Quad();
	case(2):
		return FistrElemTypeS::Quad();
	case(3):
		return FistrElemTypeS::Quad();
	case(4):
		return FistrElemTypeS::Quad();
	case(5):
		return FistrElemTypeS::Quad();
	case(6):
		return FistrElemTypeS::Quad();
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return "Unknown";
	}
}

vector<CNode*> CHexa::getFistrEdgeNode(size_t nedge)
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
		vNode.push_back(mvNode[3]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[0]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[6]);
		return vNode;

	case(7):
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[7]);
		return vNode;

	case(8):
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(9):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(10):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(11):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[6]);
		return vNode;

	case(12):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		return vNode;

	default:
		return vNode;
	}
}

string CHexa::getEdgeType(size_t nedge)
{
	if(nedge <= 12){
		return FistrElemTypeS::Beam();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CHexa::getMW3EdgeNum(size_t fistr_nedge)
{
	if(fistr_nedge <= 12 && fistr_nedge > 0){
		return fistr_nedge-1;
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);
		return 0;
	}
}

