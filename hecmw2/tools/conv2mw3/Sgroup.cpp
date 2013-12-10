/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Sgroup.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "Sgroup.h"
class CMesh;
CSgroup::CSgroup()
{
}
CSgroup::~CSgroup()
{
}
//--
//
//--
void CSgroup::addElemFaceID(size_t nElemID, size_t nFaceN)
{
	pair<size_t,size_t> nPair;

	nPair.first = nElemID;
	nPair.second= nFaceN;

	mvElemFace.push_back(nPair);
}


//--
// DATA出力準備 
//--
bool CSgroup::setup(map<size_t,CNode*> mNode, map<size_t,CElement*> mElem)
{
	vector<size_t> vNID;//節点ID

	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first; //FrontISTR 要素ID
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号
		
		CElement *pElem= mElem[nEID];
		vector<CNode*> vNode= pElem->getFistrFaceNode(nFace);//FrontISTR 面番号の節点群

		for(size_t ii=0; ii < vNode.size(); ii++){
			vNID.push_back( vNode[ii]->getID() );
		};
	};
	sort(vNID.begin(), vNID.end());
  vector<size_t>::iterator new_end = unique(vNID.begin(), vNID.end());
  vNID.erase(new_end, vNID.end());

	//出力用 Node* 配列
	for(size_t i=0; i < vNID.size(); i++){
		size_t nNID= vNID[i];
		mvNode.push_back( mNode[nNID] );
	};
	//出力用 NodeID⇒index
	for(size_t i=0; i < mvNode.size(); i++){
		CNode *pNode=mvNode[i];
		mmNID2NIX[pNode->getID()]=i;
	};
	//面の節点数
	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first; //FrontISTR 要素ID
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号

		CElement *pElem= mElem[nEID];
		mvFaceNodeNum.push_back( pElem->getFaceNodeNum(nFace) );
	};
	//面のタイプ
	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first; //FrontISTR 要素ID
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号

		CElement *pElem= mElem[nEID];
		mvFaceType.push_back( pElem->getFaceType(nFace) );
	};

	//面番号(MW3)
	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first;
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号

		CElement *pElem=mElem[nEID];

		size_t nMW3FaceN= pElem->getMW3FaceNum(nFace);
		mvMW3FaceNum.push_back(nMW3FaceN);
	};
	//面構成のBNode番号(index番号)
	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first;
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号

		CElement *pElem=mElem[nEID];

		vector<size_t> vIndex;//面構成 Node_Index
		vector<CNode*> vNode = pElem->getFistrFaceNode(nFace);//面構成Node*配列
		for(size_t ii=0; ii < vNode.size(); ii++){
			CNode* pNode=vNode[ii];
			vIndex.push_back( mmNID2NIX[pNode->getID()] );
		};
		mvvBNodeNum.push_back(vIndex);
	};

	//面構成のNodeID
	for(size_t i=0; i < mvElemFace.size(); i++){
		size_t nEID = mvElemFace[i].first;
		size_t nFace= mvElemFace[i].second;//FrontISTR 面番号

		CElement *pElem=mElem[nEID];

		vector<size_t> vNID;
		vector<CNode*> vNode = pElem->getFistrFaceNode(nFace);
		for(size_t ii=0; ii < vNode.size(); ii++){
			CNode* pNode=vNode[ii];
			vNID.push_back(pNode->getID());
		};
		mvvFaceNodeID.push_back(vNID);
	};

	return true;
}
//--
// 出力
//--
size_t CSgroup::getNumOfNode()
{
	return mvNode.size();
}
//--
CNode* CSgroup::getNode(size_t index)
{
	return mvNode[index];
}
//--
size_t CSgroup::getNumOfFace()
{
	return mvElemFace.size();
}
//--
size_t CSgroup::getElemID(size_t index)
{
	return mvElemFace[index].first;
}
//--
size_t CSgroup::getFaceN(size_t index)//面番号 入力値のまま.
{
	return mvElemFace[index].second;
}
//--
size_t CSgroup::getFaceNodeNum(size_t index)//面の節点数
{
	return mvFaceNodeNum[index];
}
//--
string CSgroup::getFaceType(size_t index)//面のタイプ
{
	return mvFaceType[index];
}
//--
size_t CSgroup::getConvMW3FaceN(size_t index)//FrontISTR⇒MW3 面番号
{
	return mvMW3FaceNum[index];
}
//--
vector<size_t> CSgroup::getBNodeN_Face(size_t index)//面を構成するBoundaryNode番号(mvNodeのindex番号)
{
	return mvvBNodeNum[index];
}
//--
vector<size_t> CSgroup::getFaceNodeID(size_t index)
{
	return mvvFaceNodeID[index];
}




