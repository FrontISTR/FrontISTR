/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   GroupPair.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "GroupPair.h"

CGroupPair::CGroupPair()
{
}
CGroupPair::~CGroupPair()
{
}
//--
//
//--
void CGroupPair::setPair(CSgroup* maSgrp, CSgroup* slSgrp)
{
	mSgrpPair.first= maSgrp;
	mSgrpPair.second= slSgrp;
}

void CGroupPair::setName(string name)
{
	msName= name;
}
string CGroupPair::getName()
{
	return msName;
}
void CGroupPair::setID(size_t id)
{
	mID=id;
}
size_t CGroupPair::getID()
{
	return mID;
}
void CGroupPair::setMeshPair(CMesh* pMMesh, CMesh* pSMesh)
{ 
	mMeshPair.first=pMMesh;
	mMeshPair.second=pSMesh;
}

//--
// DATA出力準備
//--
bool CGroupPair::setup()
{
	//マスター面群 setup
	CMesh *pMMesh= mMeshPair.first;
	map<size_t,CNode*> mMNode= pMMesh->getMapNode();
	map<size_t,CElement*> mMElem= pMMesh->getMapElement();

	mSgrpPair.first->setup(mMNode, mMElem);//---------------マスター面群 setup

	//スレーブ面群 setup
	CMesh *pSMesh= mMeshPair.second;
	map<size_t,CNode*> mSNode= pSMesh->getMapNode();
	map<size_t,CElement*> mSElem= pSMesh->getMapElement();

	mSgrpPair.second->setup(mSNode, mSElem);//--------------スレーブ面群 setup


	//----
	// MeshID,NodeID => ConIX
	//---
	mmvPartsNID2IX.resize(2);//0:マスター, 1:スレーブ
	size_t nCountIX(0);// <-----------------------ConIX

	for(size_t i=0; i < mSgrpPair.first->getNumOfNode(); i++){
		CNode* pNode= mSgrpPair.first->getNode(i);
		size_t nID= pNode->getID();
		mmvPartsNID2IX[0][nID]=nCountIX;
		nCountIX++;//-------------------------------ConIX
		mvNode.push_back(pNode);//------------------Node*

		if(i==0){ mMaxNodeID=nID; mMinNodeID=nID; }
		if(i!=0){
			if(nID > mMaxNodeID) mMaxNodeID=nID;
			if(nID < mMinNodeID) mMinNodeID=nID;
		}
	};
	for(size_t i=0; i < mSgrpPair.second->getNumOfNode(); i++){
		CNode* pNode= mSgrpPair.second->getNode(i);
		size_t nID= pNode->getID();
		mmvPartsNID2IX[1][nID]=nCountIX;
		nCountIX++;//-------------------------------ConIX
		mvNode.push_back(pNode);//------------------Node*
		if(nID > mMaxNodeID) mMaxNodeID=nID;
		if(nID < mMinNodeID) mMinNodeID=nID;
	};
	mMaxNodeIX = ( mSgrpPair.first->getNumOfNode() + mSgrpPair.second->getNumOfNode() )-1;
	mMinNodeIX = 0;

	//--
	//面構成NodeIX(ContactMesh接合節点番号)
	//--
	//マスター
	//--
	size_t nMFNum= mSgrpPair.first->getNumOfFace();
	for(size_t i=0; i < nMFNum; i++){
		vector<size_t> vNID=mSgrpPair.first->getFaceNodeID(i);
		vector<size_t> vConIX;
		for(size_t ii=0; ii < vNID.size(); ii++){
			size_t nID=vNID[ii];
			size_t nConIX= mmvPartsNID2IX[0][nID];//マスターパーツNodeID⇒ConIX
			vConIX.push_back(nConIX);
		};
		mvvMFaceNodeIX.push_back(vConIX);
	};
	//--
	//スレーブ
	//--
	size_t nSFNum= mSgrpPair.second->getNumOfFace();
	for(size_t i=0; i < nSFNum; i++){
		vector<size_t> vNID=mSgrpPair.second->getFaceNodeID(i);
		vector<size_t> vConIX;
		for(size_t ii=0; ii < vNID.size(); ii++){
			size_t nID=vNID[ii];
			size_t nConIX= mmvPartsNID2IX[1][nID];//スレーブパーツNodeID⇒ConIX
			vConIX.push_back(nConIX);
		};
		mvvSFaceNodeIX.push_back(vConIX);
	};

	mMaxMFaceN=0; mMinMFaceN=nMFNum-1;
	mMaxSFaceN=0; mMinSFaceN=nSFNum-1;

	return true;
}

size_t CGroupPair::getNumOfAllNode()
{
	size_t nAllNum= getNumOfMasterNode() + getNumOfSlaveNode();

	return nAllNum;
}
CNode* CGroupPair::getNode(size_t index)
{
	return mvNode[index];
}
size_t CGroupPair::getNumOfMasterNode()
{
	return mSgrpPair.first->getNumOfNode();
}
size_t CGroupPair::getNumOfSlaveNode()
{
	return mSgrpPair.second->getNumOfNode();
}

size_t CGroupPair::getNumOfMasterFace()
{
	return mSgrpPair.first->getNumOfFace();
}
size_t CGroupPair::getNumOfSlaveFace()
{
	return mSgrpPair.second->getNumOfFace();
}
vector<size_t> CGroupPair::getMFaceNodeIX(size_t imface)
{
	return mvvMFaceNodeIX[imface];
}
vector<size_t> CGroupPair::getSFaceNodeIX(size_t isface)
{
	return mvvSFaceNodeIX[isface];
}

size_t CGroupPair::getMasterMeshID()
{
	return mMeshPair.first->getID();
}
size_t CGroupPair::getSlaveMeshID()
{
	return mMeshPair.second->getID();
}

size_t CGroupPair::getMaxNodeIX()
{
	return mMaxNodeIX;
}
size_t CGroupPair::getMinNodeIX()
{
	return mMinNodeIX;
}

size_t CGroupPair::getMasterElemID(size_t imface)
{
	return mSgrpPair.first->getElemID(imface);
}
size_t CGroupPair::getSlaveElemID(size_t isface)
{
	return mSgrpPair.second->getElemID(isface);
}

size_t CGroupPair::getMasterFaceN_MW3(size_t imface)
{
	return mSgrpPair.first->getConvMW3FaceN(imface);
}
size_t CGroupPair::getSlaveFaceN_MW3(size_t isface)
{
	return mSgrpPair.second->getConvMW3FaceN(isface);
}

string CGroupPair::getMasterFaceType(size_t imface)
{
	return mSgrpPair.first->getFaceType(imface);
}
string CGroupPair::getSlaveFaceType(size_t isface)
{
	return mSgrpPair.second->getFaceType(isface);
}




