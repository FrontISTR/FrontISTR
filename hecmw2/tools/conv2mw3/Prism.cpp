/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Prism.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Prism.h"


CPrism::CPrism(void)
{
	mvNode.resize(6);
}
CPrism::~CPrism(void)
{
}

vector<CNode*> CPrism::getFistrFaceNode(size_t nface)
{ 
	vector<CNode*> vNode;

	switch(nface){
	case(1):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		break;
	case(2):
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		break;
	case(3):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[3]);
		break;
	case(4):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[4]);
		break;
	case(5):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[5]);
		break;
	}

	return vNode;
}
size_t CPrism::getMW3FaceNum(size_t fistr_nface)
{ 
	size_t mw_nface = fistr_nface-1;
	return mw_nface;
}

size_t CPrism::getNodeID_Fistr2MW3(size_t index)
{
	return CElement::getNodeID(index);
}

size_t CPrism::getFaceNodeNum(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return 3;
	case(2):
		return 3;
	case(3):
		return 4;
	case(4):
		return 4;
	case(5):
		return 4;
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return 0;
	}
}
string CPrism::getFaceType(size_t fistr_nface)
{
	switch(fistr_nface){
	case(1):
		return FistrElemTypeS::Triangle();
	case(2):
		return FistrElemTypeS::Triangle();
	case(3):
		return FistrElemTypeS::Quad();
	case(4):
		return FistrElemTypeS::Quad();
	case(5):
		return FistrElemTypeS::Quad();
	default:
		{CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getFaceType, arg was out of range, elem id:"+sID);}
		return "Unknown";
	}
}

vector<CNode*> CPrism::getFistrEdgeNode(size_t nedge)
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
		vNode.push_back(mvNode[3]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(5):
		vNode.push_back(mvNode[4]);
		vNode.push_back(mvNode[5]);
		return vNode;

	case(6):
		vNode.push_back(mvNode[5]);
		vNode.push_back(mvNode[3]);
		return vNode;

	case(7):
		vNode.push_back(mvNode[0]);
		vNode.push_back(mvNode[3]);
		return vNode;

	case(8):
		vNode.push_back(mvNode[1]);
		vNode.push_back(mvNode[4]);
		return vNode;

	case(9):
		vNode.push_back(mvNode[2]);
		vNode.push_back(mvNode[5]);
		return vNode;

	default:
		return vNode;
	}
}

string CPrism::getEdgeType(size_t nedge)
{
	if(nedge <= 9){
		return FistrElemTypeS::Beam();
	}else{
		CMessage *pMsg=CMessage::Instance();
		string sID=lexical_cast<string>(mID);
		pMsg->error("getEdgeType, arg was out of range, elem id:"+sID);

		return FistrElemTypeS::Unknown();
	}
}

size_t CPrism::getMW3EdgeNum(size_t fistr_nedge)
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




