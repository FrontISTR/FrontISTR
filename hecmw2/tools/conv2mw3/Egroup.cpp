/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Egroup.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#include "Egroup.h"

CEgroup::CEgroup()
{
}
CEgroup::~CEgroup()
{
}


//--
// 
//--
void CEgroup::addElement(CElement* pElem)
{
	mvEID.push_back(pElem->getID());
	
	mvTempElem.push_back(pElem);
	mmEID2IX[ pElem->getID() ] = mvTempElem.size()-1;
}

//--
// DATA出力準備
//--
bool CEgroup::setup()
{
	// 要素セットアップ
	sort(mvEID.begin(), mvEID.end());
	vector<size_t>::iterator new_eend = unique(mvEID.begin(), mvEID.end());
	mvEID.erase(new_eend, mvEID.end());

	for(size_t i=0; i < mvEID.size(); i++){
		size_t nID= mvEID[i];
		size_t index= mmEID2IX[nID];
		CElement *pElem=mvTempElem[index];

		mvElement.push_back(pElem);
	};


	// Node番号セットアップ
	for(size_t i=0; i < mvElement.size(); i++){
		CElement *pElem=mvElement[i];

		size_t nNNum= pElem->getNumOfNode();
		for(size_t ii=0; ii < nNNum; ii++){
			size_t nID = pElem->getNodeID(ii);
			mvNID.push_back(nID);
		};
	};

	sort(mvNID.begin(), mvNID.end());
  vector<size_t>::iterator new_nend = unique(mvNID.begin(), mvNID.end());
  mvNID.erase(new_nend, mvNID.end());

	for(size_t i=0; i < mvNID.size(); i++){
		size_t nID= mvNID[i];
		mmNID2IX[nID] = i;
	}

	return true;
}

size_t CEgroup::getNumOfNode()
{
	return mvNID.size();
}
size_t CEgroup::getNodeID(size_t index)
{
	return mvNID[index];
}
size_t CEgroup::getBNodeIndex(size_t nid)
{
	return mmNID2IX[nid];
}




