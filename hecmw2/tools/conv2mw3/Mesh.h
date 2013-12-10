/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Mesh.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef E6805807_B095_41e3_B709_86DE26F977CC
#define E6805807_B095_41e3_B709_86DE26F977CC

#include <vector>
#include <map>
#include <algorithm>
using namespace std;

#include <boost/lexical_cast.hpp>
using namespace boost;

#include "Node.h"
#include "Element.h"

#include "Beam.h"
#include "Beam2.h"
#include "Hexa.h"
#include "Hexa2.h"
#include "IFaceQuad.h"
#include "IFaceQuad2.h"
#include "Prism.h"
#include "Prism2.h"
#include "Quad.h"
#include "Quad2.h"
#include "QuadShell.h"
#include "QuadShell2.h"
#include "Tetra.h"
#include "Tetra2.h"
#include "Triangle.h"
#include "Triangle2.h"
#include "TriShell.h"
#include "TriShell2.h"

#include "Group.h"
#include "Ngroup.h"
#include "Sgroup.h"
#include "Egroup.h"
#include "Lgroup.h"

#include "Message.h"

class CMesh{
public:
	CMesh();
	virtual ~CMesh();

private:
	size_t mID;

	vector<CNode*>       mvNode;
	vector<CElement*> mvElement;

	map<size_t, CNode*> mmNode;
	map<size_t, CElement*> mmElement;

	vector<size_t>   mvNID,  mvEID;//-- IDî‘çÜóÒ
	size_t mMaxNID, mMinNID, mMaxEID, mMinEID;
	
	string msName;//-- parts name

	// Grp
	vector<CNgroup*> mvNgrp;
	vector<CSgroup*> mvSgrp;
	vector<CEgroup*> mvEgrp;
	vector<CLgroup*> mvLgrp;

	map<string, CNgroup*> mmNgrp;
	map<string, CSgroup*> mmSgrp;
	map<string, CEgroup*> mmEgrp;
	map<string, CLgroup*> mmLgrp;

	vector<size_t> mvNgrpID, mvSgrpID, mvEgrpID, mvLgrpID;//-- IDî‘çÜóÒ

	map<string, size_t> mNgrpName2IX;//GroupNameÅÀNgrpÇÃIndex
	map<string, size_t> mSgrpName2IX;//GroupNameÅÀSgrpÇÃIndex
	map<string, size_t> mEgrpName2IX;//GroupNameÅÀEgrpÇÃIndex
	map<string, size_t> mLgrpName2IX;//GropuNameÅÀLgrpÇÃIndex


	//åüçıï∂éöóÒ
	string msNgrpNameS;
	string msSgrpNameS;
	string msEgrpNameS;
	string msLgrpNameS;

private:
	bool existID(vector<size_t>& vec, size_t id);

public:
	void setPartsName(string sname);
	void setID(size_t id);
	size_t getID(){ return mID;}
	string getName(){ return msName;}

	void addNode(CNode* pNode);
	void addElement(CElement* pElement);

	CNode* getNode(size_t index);
	CNode* getNode_id(size_t id);
	CElement* getElement(size_t index);
	CElement* getElement_id(size_t id);

	size_t getNumOfNode(){ return mvNode.size();}
	size_t getNumOfElement(){ return mvElement.size();}

	size_t getMaxNodeID(){ return mMaxNID;}
	size_t getMinNodeID(){ return mMinNID;}
	size_t getMaxElementID(){ return mMaxEID;}
	size_t getMinElementID(){ return mMinEID;}

	//--
	// id åüçı
	//--
	bool existNodeID(size_t id);
	bool existElementID(size_t id);


	// GroupñºëOåüçı
	bool existNgrpName(string sNgrpName);
	bool existSgrpName(string sSgrpName);
	bool existEgrpName(string sEgrpName);
	bool existLgrpName(string sLgrpName);

	// Group
	void addNgroup(CNgroup* pNgrp);
	void addSgroup(CSgroup* pSgrp);
	void addEgroup(CEgroup* pEgrp);
	void addLgroup(CLgroup* pLgrp);

	size_t getNumOfNgrp();
	CNgroup* getNgrp(size_t index);
	CNgroup* getNgrp(string sNgrpName);

	size_t getNumOfSgrp();
	CSgroup* getSgrp(size_t index);
	CSgroup* getSgrp(string sSgrpName);

	size_t getNumOfEgrp();
	CEgroup* getEgrp(size_t index);
	CEgroup* getEgrp(string sEgrpName);

	size_t getNumOfLgrp();
	CLgroup* getLgrp(size_t index);
	CLgroup* getLgrp(string sLgrpName);

	//--
	// DATAèoóÕèÄîı
	//--
	bool setup();

	// èoóÕêßå‰
	size_t getNgrpIndex(string sNgrpName);
	size_t getSgrpIndex(string sSgrpName);
	size_t getEgrpIndex(string sEgrpName);
	size_t getLgrpIndex(string sLgrpName);

	// for ê⁄çáñ setup
	map<size_t, CNode*> getMapNode(){ return mmNode;}
	map<size_t, CElement*> getMapElement(){ return mmElement;}

};
#endif //include_guard



















