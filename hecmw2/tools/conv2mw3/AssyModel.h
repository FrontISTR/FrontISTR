/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   AssyModel.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/

#ifndef DE4C3D2B_A046_4ecc_B8FC_C0C83F11BECD
#define DE4C3D2B_A046_4ecc_B8FC_C0C83F11BECD

#include "Mesh.h"
#include <vector>
#include <map>
#include <algorithm>
#include <string>
using namespace std;

#include <boost/lexical_cast.hpp>
using namespace boost;

#include "GroupPair.h"
#include "MPCPair.h" //接合ペア
#include "ConPair.h" //接触ペア

#include "Message.h"

class CAssyModel{
public:
	CAssyModel();
	virtual ~CAssyModel();

private:
	// Mesh:パーツ
	vector<CMesh*> mvMesh;
	map<string,CMesh*> mmMesh;
	vector<size_t> mvMID;//MeshのID番号管理
	size_t mMaxMID,mMinMID;

	// GroupPair
	vector<CGroupPair*> mvMPCPair;//接合
	map<string, CGroupPair*> mmMPCPair;
	vector<size_t> mvMPCID;//MPCペアのID番号管理
	size_t mMaxMPCID, mMinMPCID;

	vector<CGroupPair*> mvConPair;//接触
	map<string, CGroupPair*> mmConPair;
	vector<size_t> mvConID;//ContactペアのID番号管理
	size_t mMaxConID, mMinConID;


	//検索文字列(パーツ名, 接合面名)
	string msPartsNameS;
	string msMPCPairNameS;
	string msConPairNameS;

public:
	//--
	// Mesh(Egrp,Sgrp,Ngrp)
	//--
	void addMesh(CMesh* pMesh);

	size_t getNumOfMesh();
	CMesh* getMesh(size_t index);
	CMesh* getMesh(string sPartsName);

	//--
	// 名前検索
	//--
	bool existPartsName(string sPartsName);
	bool existMPCPairName(string sPairName);
	bool existConPairName(string sPairName);

	//--
	// GroupPair
	//--
	void addMPCPair(CGroupPair* pPair);
	void addConPair(CGroupPair* pPair);

	size_t getNumOfMPCPair();
	CGroupPair* getMPCPair(size_t index);
	CGroupPair* getMPCPair(string sMPCName);

	size_t getNumOfConPair();
	CGroupPair* getConPair(size_t index);
	CGroupPair* getConPair(string sConName);

public:
	//--
	// DATA出力準備
	//--
	bool setup();

	bool existMesh();
	size_t getMaxIDinMesh();
	size_t getMinIDinMesh();

	bool existMPCPair();
	size_t getMaxIDinMPCPair();
	size_t getMinIDinMPCPair();

	bool existConPair();
	size_t getMaxIDinConPair();
	size_t getMinIDinConPair();
};

#endif //include_guard


