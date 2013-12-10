/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Tetra2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Tetra2.h"


CTetra2::CTetra2(void)
{
	mvNode.resize(10);
}
CTetra2::~CTetra2(void)
{
}

vector<CNode*> CTetra2::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		break;
	case(2):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[8]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		break;
	case(3):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[8]);
		break;
	case(4):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		break;
	}

	return vNode;
}
size_t CTetra2::getMW3FaceNum(size_t fistr_nface)
{ 
	size_t mw_nface = fistr_nface-1;
	return mw_nface;
}

size_t CTetra2::getNodeID_Fistr2MW3(size_t index)
{
	CMessage *pMsg=CMessage::Instance();

	size_t up_index;

	switch(index){
	case(0):up_index= 0; break;
	case(1):up_index= 1; break;
	case(2):up_index= 2; break;
	case(3):up_index= 3; break;
	case(4):up_index= 5; break;
	case(5):up_index= 6; break;
	case(6):up_index= 4; break;
	case(7):up_index= 7; break;
	case(8):up_index= 8; break;
	case(9):up_index= 9; break;
	default:
		{string sIX = lexical_cast<string>(index);
		pMsg->error(sIX+" was out of range, Tetra2::getNodeID_Fistr2MW3");
		pMsg->info("index set to 9");}
		up_index= 9;
		break;
	}

	return mvNode[up_index]->getID();
}

size_t CTetra2::getFaceNodeNum(size_t fistr_nface)
{
	if(fistr_nface > 4){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return 0;
	}

	return 6;
}
string CTetra2::getFaceType(size_t fistr_nface)
{
	if(fistr_nface > 4){
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);
		return "Unknown";
	}

	return FistrElemTypeS::Triangle2();
}

vector<CNode*> CTetra2::getFistrEdgeNode(size_t nedge)
{
	vector<CNode*> vNode;
	//--
	//FrontISTR‚Í"1"‚©‚çŽn‚Ü‚é
	//--
	switch(nedge){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[6]);
		return vNode;

	case(2):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[7]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[8]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[9]);
		return vNode;

	default:
		return vNode;
	}
}

string CTetra2::getEdgeType(size_t nedge)
{
	if(nedge <= 6){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CTetra2::getMW3EdgeNum(size_t fistr_nedge)
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


