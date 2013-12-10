/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Prism2.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Prism2.h"


CPrism2::CPrism2(void)
{
	mvNode.resize(15);
}

CPrism2::~CPrism2(void)
{
}

vector<CNode*> CPrism2::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[8]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[7]);
		break;
	case(2):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[11]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[10]);
		break;
	case(3):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[8]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[13]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[11]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[12]);
		break;
	case(4):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[6]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[14]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[9]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[13]);
		break;
	case(5):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[7]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[12]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[10]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[14]);
		break;
	}

	return vNode;
}
size_t CPrism2::getMW3FaceNum(size_t fistr_nface)
{ 
	size_t mw_nface = fistr_nface-1;
	return mw_nface;
}

size_t CPrism2::getNodeID_Fistr2MW3(size_t index)
{
	CMessage *pMsg=CMessage::Instance();

	size_t up_index;

	switch(index){
	case(0):up_index= 0; break;
	case(1):up_index= 1; break;
	case(2):up_index= 2; break;
	case(3):up_index= 3; break;
	case(4):up_index= 4; break;
	case(5):up_index= 5; break;

	case(6):up_index= 8; break;
	case(7):up_index= 7; break;
	case(8):up_index= 6; break;

	case(9): up_index= 13; break;
	case(10):up_index= 14; break;
	case(11):up_index= 12; break;

	case(12):up_index=  9; break;
	case(13):up_index= 10; break;
	case(14):up_index= 11; break;

	default:
		string sIX = lexical_cast<string>(index);
		pMsg->error(sIX+" was out of range, Prism2::getNodeID_Fistr2MW3");
		pMsg->info("index set to 14");

		up_index= 14;
		break;
	}

	return mvNode[up_index]->getID();
}

size_t CPrism2::getFaceNodeNum(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return 6;
	case(2):
		return 6;
	case(3):
		return 8;
	case(4):
		return 8;
	case(5):
		return 8;
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return 0;
	}
}
string CPrism2::getFaceType(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return FistrElemTypeS::Triangle2();
	case(2):
		return FistrElemTypeS::Triangle2();
	case(3):
		return FistrElemTypeS::Quad2();
	case(4):
		return FistrElemTypeS::Quad2();
	case(5):
		return FistrElemTypeS::Quad2();
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return "Unknown";
	}
}

vector<CNode*> CPrism2::getFistrEdgeNode(size_t nedge)
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
		vNode.push_back(mvNode[6]);
		return vNode;

	case(3):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[7]);
		return vNode;

	case(4):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[11]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[9]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[10]);
		return vNode;

	case(7):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[12]);
		return vNode;

	case(8):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[13]);
		return vNode;

	case(9):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[14]);
		return vNode;

	default:
		return vNode;
	}
}

string CPrism2::getEdgeType(size_t nedge)
{
	if(nedge <= 9){
		return FistrElemTypeS::Beam2();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CPrism2::getMW3EdgeNum(size_t fistr_nedge)
{
	switch(fistr_nedge){
	case(1):
		return 0;
	case(2):
		return 2;
	case(3):
		return 1;

	case(4):
		return 6;
	case(5):
		return 7;
	case(6):
		return 8;

	case(7):
		return 3;
	case(8):
		return 4;
	case(9):
		return 5;

	default:
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getMW3EdgeNum, arg was out of range, elem id:"+sID);
		return 0;
	}
}


