/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Hexa2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Hexa2.h"


CHexa2::CHexa2(void)
{
	mvNode.resize(20);
}

CHexa2::~CHexa2(void)
{
}

vector<CNode*> CHexa2::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[8]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[10]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[11]);
		break;
	case(2):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[12]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[13]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[14]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[15]);
		break;
	case(3):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[8]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[17]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[12]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[16]);
		break;
	case(4):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[18]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[13]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[17]);
		break;
	case(5):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[10]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[19]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[14]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[18]);
		break;
	case(6):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[11]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[16]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[15]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[19]);
		break;
	}

	return vNode;
}
size_t CHexa2::getMW3FaceNum(size_t fistr_nface)
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

size_t CHexa2::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}


size_t CHexa2::getFaceNodeNum(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return 8;
	case(2):
		return 8;
	case(3):
		return 8;
	case(4):
		return 8;
	case(5):
		return 8;
	case(6):
		return 8;
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return 0;
	}
}
string CHexa2::getFaceType(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return FistrElemTypeS::Quad2();
	case(2):
		return FistrElemTypeS::Quad2();
	case(3):
		return FistrElemTypeS::Quad2();
	case(4):
		return FistrElemTypeS::Quad2();
	case(5):
		return FistrElemTypeS::Quad2();
	case(6):
		return FistrElemTypeS::Quad2();
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return "Unknown";
	}
}

vector<CNode*> CHexa2::getFistrEdgeNode(size_t nedge)
{
	vector<CNode*> vNode;
	//--
	//FrontISTR‚Í"1"‚©‚çŽn‚Ü‚é
	//--
	switch(nedge){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[8]);
		return vNode;

	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[9]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[10]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[11]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[12]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[13]);
		return vNode;

	case(7):
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[14]);
		return vNode;

	case(8):
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[15]);
		return vNode;

	case(9):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[16]);
		return vNode;

	case(10):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[17]);
		return vNode;

	case(11):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[18]);
		return vNode;

	case(12):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[19]);
		return vNode;

	default:
		return vNode;
	}
}

string CHexa2::getEdgeType(size_t nedge)
{
	if(nedge <= 12){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CHexa2::getMW3EdgeNum(size_t fistr_nedge)
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





