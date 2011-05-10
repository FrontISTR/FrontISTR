//  ContactElement.h
//  (ContactElement:Face)-> Face
//
//				2009.01.08
//				2009.01.08
//				k.Takeda
#ifndef CONTACT_ELEMENT_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC
#define CONTACT_ELEMENT_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC

#include "CommonStd.h"
#include "Node.h"

namespace pmw{
class CContactElement{
public:
	CContactElement();

virtual ~CContactElement();


protected:
	uint mnAssyID;//

	// Face
	vector<CNode*> mvNode;

public:
	// Face.push_back
	void addNode(CNode* pNode){ mvNode.push_back(pNode);}

	// Face":=" setNode(vNode)
	void operator=(vector<CNode*> vNode){ mvNode = vNode;}


	// Face
	void resizeNode(const uint& size){  mvNode.resize(size);}
	void setNode(CNode* pNode, const uint& index){ mvNode[index] = pNode;}
	void setNode(vector<CNode*>& vNode){ mvNode = vNode;}// same "="


	//
	uint   getNumOfNode() { return mvNode.size();}
	vector<CNode*>&  getNodes(){ return mvNode;}
	CNode*  getNode(const uint& i){ return mvNode[i];}


	//
	void  setAssyID(const uint& id){ mnAssyID = id;}
	uint& getAssyID(){ return mnAssyID;}
};
}
#endif

