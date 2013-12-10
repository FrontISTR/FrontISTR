/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Mesh.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#include "Mesh.h"

CMesh::CMesh()
{
	msName="";//parts name

	msNgrpNameS="";
	msSgrpNameS="";
	msEgrpNameS="";
}
CMesh::~CMesh()
{
	vector<CNode*>::iterator itn;
	for(itn=mvNode.begin(); itn != mvNode.end(); ++itn)  delete *itn;

	vector<CElement*>::iterator ite;
	for(ite=mvElement.begin(); ite != mvElement.end(); ++ite)  delete *ite;
	
	vector<CNgroup*>::iterator itng;
	for(itng=mvNgrp.begin(); itng != mvNgrp.end(); ++itng) delete *itng;

	vector<CSgroup*>::iterator itsg;
	for(itsg=mvSgrp.begin(); itsg != mvSgrp.end(); ++itsg) delete *itsg;

	vector<CEgroup*>::iterator iteg;
	for(iteg=mvEgrp.begin(); iteg != mvEgrp.end(); ++iteg) delete *iteg;
}
//--
//
//--
void CMesh::setPartsName(string sname)
{
	msName= sname;
}
void CMesh::setID(size_t id)
{
	mID= id;
}
//--
//
//--
void CMesh::addNode(CNode* pNode)
{
	mvNode.push_back(pNode);
	mmNode[pNode->getID()]= pNode;

	mvNID.push_back(pNode->getID());//-- ID列
}
void CMesh::addElement(CElement* pElement)
{
	mvElement.push_back(pElement);
	mmElement[pElement->getID()]=pElement;

	mvEID.push_back(pElement->getID());//-- ID列
}
//--
//
//--
CNode* CMesh::getNode(size_t index)
{
	return mvNode[index];
}
CNode* CMesh::getNode_id(size_t id)
{
	return mmNode[id];
}
CElement* CMesh::getElement(size_t index)
{
	return mvElement[index];
}
CElement* CMesh::getElement_id(size_t id)
{
	return mmElement[id];
}
//--
// id 検索
//--
bool CMesh::existID(vector<size_t>& vec, size_t id)
{
	// vec初期状態
	if(vec.empty()) return false;
	
	size_t low=0, high=vec.size()-1, ix;
	// 検索
	while(low <= high){

		ix=(low+high)/2;

		if(id == vec[ix]){
			return true;//------------ true
		}else if(id < vec[ix]){
			high= ix-1;
		}else{
			low = ix+1;
		}
	};//while end

	return false;
}
bool CMesh::existNodeID(size_t id)
{
	return existID(mvNID, id);
}
bool CMesh::existElementID(size_t id)
{
	return existID(mvEID, id);
}
//--
// Group名検索
//--
bool CMesh::existNgrpName(string sNgrpName)
{
	string::size_type pos;

	pos= msNgrpNameS.find(sNgrpName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
bool CMesh::existSgrpName(string sSgrpName)
{
	string::size_type pos;

	pos= msSgrpNameS.find(sSgrpName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
bool CMesh::existEgrpName(string sEgrpName)
{
	string::size_type pos;

	pos= msEgrpNameS.find(sEgrpName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
bool CMesh::existLgrpName(string sLgrpName)
{
	string::size_type pos;

	pos= msLgrpNameS.find(sLgrpName);
	if(pos!=string::npos){
		return true;
	}else{
		return false;
	}
}
//--
// add Group
//--
void CMesh::addNgroup(CNgroup* pNgrp)
{
	string sGrpName= pNgrp->getGroupName();

	mvNgrp.push_back(pNgrp);
	mmNgrp[sGrpName]= pNgrp;

	msNgrpNameS += sGrpName + ' ';
}
void CMesh::addSgroup(CSgroup* pSgrp)
{
	string sSgrpName= pSgrp->getGroupName();

	mvSgrp.push_back(pSgrp);
	mmSgrp[sSgrpName]= pSgrp;

	msSgrpNameS += sSgrpName + ' ';
}
void CMesh::addEgroup(CEgroup* pEgrp)
{
	string sEgrpName= pEgrp->getGroupName();

	mvEgrp.push_back(pEgrp);
	mmEgrp[sEgrpName]= pEgrp;

	msEgrpNameS += sEgrpName + ' ';
}
void CMesh::addLgroup(CLgroup* pLgrp)
{
	string sLgrpName= pLgrp->getGroupName();

	mvLgrp.push_back(pLgrp);
	mmLgrp[sLgrpName]= pLgrp;

	msLgrpNameS += sLgrpName + ' ';
}

//--
// get Group
//--
size_t CMesh::getNumOfNgrp()
{
	return mvNgrp.size();
}
CNgroup* CMesh::getNgrp(size_t index)
{
	return mvNgrp[index];
}
CNgroup* CMesh::getNgrp(string sNgrpName)
{
	return mmNgrp[sNgrpName];
}
//----
size_t CMesh::getNumOfSgrp()
{
	return mvSgrp.size();
}
CSgroup* CMesh::getSgrp(size_t index)
{
	return mvSgrp[index];
}
CSgroup* CMesh::getSgrp(string sSgrpName)
{
	return mmSgrp[sSgrpName];
}
//----
size_t CMesh::getNumOfEgrp()
{
	return mvEgrp.size();
}
CEgroup* CMesh::getEgrp(size_t index)
{
	return mvEgrp[index];
}
CEgroup* CMesh::getEgrp(string sEgrpName)
{
	return mmEgrp[sEgrpName];
}
//----
size_t CMesh::getNumOfLgrp()
{
	return mvLgrp.size();
}
CLgroup* CMesh::getLgrp(size_t index)
{
	return mvLgrp[index];
}
CLgroup* CMesh::getLgrp(string sLgrpName)
{
	return mmLgrp[sLgrpName];
}

//--
// DATA出力準備
//--
bool CMesh::setup()
{
	// Node ID番号:最大、最小
	sort(mvNID.begin(),mvNID.end());

	size_t nsize=mvNID.size();
	if(nsize > 0){
		mMaxNID=mvNID[nsize-1];
		mMinNID=mvNID[0];
	}
	// Element ID番号:最大、最小
	sort(mvEID.begin(),mvEID.end());

	nsize=mvEID.size();
	if(nsize > 0){
		mMaxEID=mvEID[nsize-1];
		mMinEID=mvEID[0];
	}

	size_t id;
	// Ngrp ID番号
	for(size_t i=0; i < mvNgrp.size(); i++){
		CNgroup *pNgrp=mvNgrp[i];
		id= pNgrp->getID();
		mvNgrpID.push_back(id);
	};
	sort(mvNgrpID.begin(), mvNgrpID.end());

	// Sgrp ID番号
	for(size_t i=0; i < mvSgrp.size(); i++){
		CSgroup *pSgrp=mvSgrp[i];
		id= pSgrp->getID();
		mvSgrpID.push_back(id);
	};
	sort(mvSgrpID.begin(),mvSgrpID.end());

	// Egrp ID番号
	for(size_t i=0; i < mvEgrp.size(); i++){
		CEgroup *pEgrp=mvEgrp[i];
		id= pEgrp->getID();
		mvEgrpID.push_back(id);
	};
	sort(mvEgrpID.begin(),mvEgrpID.end());

	// Lgrp ID番号
	for(size_t i=0; i < mvLgrp.size(); i++){
		CLgroup *pLgrp=mvLgrp[i];
		id= pLgrp->getID();
		mvLgrpID.push_back(id);
	};
	sort( mvLgrpID.begin(), mvLgrpID.end() );



	// Ngrp 名前⇒index
	for(size_t i=0; i < mvNgrp.size(); i++){
		CNgroup *pNgrp=mvNgrp[i];
		mNgrpName2IX[pNgrp->getGroupName()]=i;
	};
	// Sgrp 名前⇒index
	for(size_t i=0; i < mvSgrp.size(); i++){
		CSgroup *pSgrp=mvSgrp[i];
		mSgrpName2IX[pSgrp->getGroupName()]=i;
	};
	// Egrp 名前⇒index
	for(size_t i=0; i < mvEgrp.size(); i++){
		CEgroup *pEgrp=mvEgrp[i];
		mEgrpName2IX[pEgrp->getGroupName()]=i;
	};
	// Lgrp 名前⇒index
	for(size_t i=0; i < mvLgrp.size(); i++){
		CLgroup *pLgrp=mvLgrp[i];
		mLgrpName2IX[pLgrp->getGroupName()]=i;
	};


	//--
	// Sgrp(面群) setup
	//--
	for(size_t i=0; i < mvSgrp.size(); i++){
		CSgroup *pSgrp=mvSgrp[i];
		pSgrp->setup(mmNode, mmElement);
	};
	//--
	// Lgrp(辺群) setup
	//--
	for(size_t i=0; i < mvLgrp.size(); i++){
		CLgroup* pLgrp=mvLgrp[i];
		pLgrp->setup(mmNode,mmElement);
	};
	//--
	// Egrp(要素群) setup
	//--
	for(size_t i=0; i < mvEgrp.size(); i++){
		CEgroup *pEgrp=mvEgrp[i];
		pEgrp->setup();
	};

	return true;
}

//--
// 出力制御
//--
size_t CMesh::getNgrpIndex(string sNgrpName)
{
	return mNgrpName2IX[sNgrpName];
}
size_t CMesh::getSgrpIndex(string sSgrpName)
{
	return mSgrpName2IX[sSgrpName];
}
size_t CMesh::getEgrpIndex(string sEgrpName)
{
	return mEgrpName2IX[sEgrpName];
}
size_t CMesh::getLgrpIndex(string sLgrpName)
{
	return mLgrpName2IX[sLgrpName];
}



