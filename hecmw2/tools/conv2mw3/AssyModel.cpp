/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   AssyModel.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "AssyModel.h"

CAssyModel::CAssyModel()
{
	msPartsNameS="";
}
CAssyModel::~CAssyModel()
{
	vector<CMesh*>::iterator it;
	for(it=mvMesh.begin(); it!=mvMesh.end(); it++) delete *it;

	vector<CGroupPair*>::iterator itmpc;
	for(itmpc=mvMPCPair.begin(); itmpc!=mvMPCPair.end(); itmpc++) delete *itmpc;

	vector<CGroupPair*>::iterator itcon;
	for(itcon=mvConPair.begin(); itcon!=mvConPair.end(); itcon++) delete *itcon;
}
//--
// Mesh
//--
void CAssyModel::addMesh(CMesh* pMesh)
{
	string sPartsName= pMesh->getName();

	mvMesh.push_back(pMesh); //vec
	mmMesh[sPartsName]=pMesh;//map

	msPartsNameS += sPartsName + ' ';
}
size_t CAssyModel::getNumOfMesh()
{
	return mvMesh.size();
}
CMesh* CAssyModel::getMesh(size_t index)
{
	return mvMesh[index];
}
CMesh* CAssyModel::getMesh(string sPartsName)
{
	return mmMesh[sPartsName];
}
//--
// 名前 検索
//--
bool CAssyModel::existPartsName(string sPartsName)
{
	string::size_type pos;

	pos= msPartsNameS.find(sPartsName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}

bool CAssyModel::existMPCPairName(string sPairName)
{
	string::size_type pos;

	pos= msMPCPairNameS.find(sPairName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
bool CAssyModel::existConPairName(string sPairName)
{
	string::size_type pos;

	pos= msConPairNameS.find(sPairName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
//--
// add GroupPair : MPC & Contact
//--
void CAssyModel::addMPCPair(CGroupPair* pPair)
{
	string sMPCName= pPair->getName();

	mvMPCPair.push_back(pPair);
	mmMPCPair[sMPCName]= pPair;

	msMPCPairNameS += sMPCName + ' ';
}
void CAssyModel::addConPair(CGroupPair* pPair)
{
	string sConName= pPair->getName();

	mvConPair.push_back(pPair);
	mmConPair[sConName]= pPair;

	msConPairNameS += sConName + ' ';
}
//--
// get GroupPair : MPC & Contact
//--
size_t CAssyModel::getNumOfMPCPair()
{
	return mvMPCPair.size();
}
CGroupPair* CAssyModel::getMPCPair(size_t index)
{
	return mvMPCPair[index];
}
CGroupPair* CAssyModel::getMPCPair(string sMPCName)
{
	return mmMPCPair[sMPCName];
}
size_t CAssyModel::getNumOfConPair()
{
	return mvConPair.size();
}
CGroupPair* CAssyModel::getConPair(size_t index)
{
	return mvConPair[index];
}
CGroupPair* CAssyModel::getConPair(string sConName)
{
	return mmConPair[sConName];
}

//--
// DATA出力準備
//--
bool CAssyModel::setup()
{
	CMessage *pMsg=CMessage::Instance();

	size_t id;
	//--
	//1.Mesh ID番号 & Mesh::setup()
	//--
	for(size_t i=0; i < mvMesh.size(); i++){
		CMesh *pMesh=mvMesh[i];

		pMesh->setup();//-------- MeshのDATA出力準備, Sgrp->setup

		id= pMesh->getID();
		mvMID.push_back(id);
	};
	sort(mvMID.begin(), mvMID.end());

	size_t msize=mvMID.size();
	if(msize > 0){
		mMaxMID=mvMID[msize-1];
		mMinMID=mvMID[0];
	}
	
	//--
	//2.MPCPair ID番号
	//--
	for(size_t i=0; i < mvMPCPair.size(); i++){
		CGroupPair *pPair=mvMPCPair[i];

		pPair->setup();//----------------GroupPairのDATA出力準備(MPC)

		id= pPair->getID();
		mvMPCID.push_back(id);
	};
	sort(mvMPCID.begin(),mvMPCID.end());

	msize=mvMPCID.size();
	if(msize > 0){
		mMaxMPCID=mvMPCID[msize-1];
		mMinMPCID=mvMPCID[0];
	}

	//--
	//3.ConPair ID番号
	//--
	for(size_t i=0; i < mvConPair.size(); i++){
		CGroupPair *pPair=mvConPair[i];

		pPair->setup();//----------------GroupPairのDATA出力準備(Contact)

		id= pPair->getID();
		mvConID.push_back(id);
	};
	sort(mvConID.begin(),mvConID.end());

	msize=mvConID.size();
	if(msize > 0){
		mMaxConID=mvConID[msize-1];
		mMinConID=mvConID[0];
	}

	pMsg->info("terminated data setup.\n");//------ メッセージ

	return true;
}
bool CAssyModel::existMesh()
{
	if(mvMesh.empty()){
		return false;
	}else{
		return true;
	}
}
size_t CAssyModel::getMaxIDinMesh()
{
	return mMaxMID;
}
size_t CAssyModel::getMinIDinMesh()
{
  return mMinMID;
}
bool CAssyModel::existMPCPair()
{
	if(mvMPCPair.empty()){
		return false;
	}else{
		return true;
	}
}
size_t CAssyModel::getMaxIDinMPCPair()
{
	return mMaxMPCID;
}
size_t CAssyModel::getMinIDinMPCPair()
{
	return mMinMPCID;
}
bool CAssyModel::existConPair()
{
	if(mvConPair.empty()){
		return false;
	}else{
		return true;
	}
}
size_t CAssyModel::getMaxIDinConPair()
{
	return mMaxConID;
}
size_t CAssyModel::getMinIDinConPair()
{
	return mMinConID;
}




