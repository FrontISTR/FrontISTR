/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   GroupPair.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include <utility>
#include "Sgroup.h"

#ifndef B8CB4CD9_C7D0_4ec1_8D86_AE207FD90EC0
#define B8CB4CD9_C7D0_4ec1_8D86_AE207FD90EC0

#include "Mesh.h"

class CGroupPair{
public:
	CGroupPair();
	virtual ~CGroupPair();

protected:
	pair<CSgroup*, CSgroup*> mSgrpPair;//Master-Slave
	pair<CMesh*, CMesh*>     mMeshPair;//Master-Slave

	string msTypeName;//-- Type Name, MPCPair || ConPair
	string msName;//------ グループ･ペア名

	size_t mID;

	vector<map<size_t,size_t> > mmvPartsNID2IX;//MeshID,NodeID => ConIX
	vector<CNode*> mvNode;//総Node*
	size_t mMaxNodeID,mMinNodeID;
	size_t mMaxNodeIX,mMinNodeIX;

	vector<vector<size_t> > mvvMFaceNodeIX;//マスター面 構成NodeIX
	vector<vector<size_t> > mvvSFaceNodeIX;//スレーブ面 構成NodeIX

	size_t mMaxMFaceN,mMinMFaceN;
	size_t mMaxSFaceN,mMinSFaceN;

public:
	void setPair(CSgroup* maSgrp, CSgroup* slSgrp);

	string getType(){ return msTypeName;}

	void setName(string name);
	string getName();

	void setID(size_t id);
	size_t getID();

	void setMeshPair(CMesh* pMMesh, CMesh* pSMesh);


	//--
	// DATA出力準備
	//--
	bool setup();

	size_t getNumOfAllNode();
	CNode* getNode(size_t index);
	size_t getMaxNodeID(){ return mMaxNodeID;}
	size_t getMinNodeID(){ return mMinNodeID;}
	size_t getMaxNodeIX();
	size_t getMinNodeIX();

	size_t getNumOfMasterNode();
	size_t getNumOfSlaveNode();

	size_t getNumOfMasterFace();
	size_t getNumOfSlaveFace();
	size_t getMaxMFaceIX(){ return mMaxMFaceN;}
	size_t getMinMFaceIX(){ return mMinMFaceN;}
	size_t getMaxSFaceIX(){ return mMaxSFaceN;}
	size_t getMinSFaceIX(){ return mMinSFaceN;}

	vector<size_t> getMFaceNodeIX(size_t imface);
	vector<size_t> getSFaceNodeIX(size_t isface);

	size_t getMasterMeshID();
	size_t getSlaveMeshID();

	size_t getMasterElemID(size_t imface);
	size_t getSlaveElemID(size_t isface);
	
	size_t getMasterFaceN_MW3(size_t imface);
	size_t getSlaveFaceN_MW3(size_t isface);

	string getMasterFaceType(size_t imface);
	string getSlaveFaceType(size_t isface);
};
#endif //include_guard

