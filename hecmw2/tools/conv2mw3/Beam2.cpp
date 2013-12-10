/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Beam2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Beam2.h"


CBeam2::CBeam2(void)
{
	mvNode.resize(3);
}
CBeam2::~CBeam2(void)
{
}

size_t CBeam2::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}

size_t CBeam2::getFaceNodeNum(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceNodeNum, arg was out of range, elem id:"+sID);
		return 0;
	}

	return 3;
}
string CBeam2::getFaceType(size_t fistr_nface)
{
	if(fistr_nface != 1){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::Beam2();
}

vector<CNode*> CBeam2::getFistrEdgeNode(size_t nedge)
{
	return mvNode;
}

string CBeam2::getEdgeType(size_t nedge)
{
	if(nedge == 1){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CBeam2::getMW3EdgeNum(size_t fistr_nedge)
{
	if(fistr_nedge == 1){
		return 0;
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);

		return 0;
	}
}

