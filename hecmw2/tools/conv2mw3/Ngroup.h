/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Ngroup.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef F2F0DD80_D3D5_4ac3_BE4C_BF29F2E7AAE7
#define F2F0DD80_D3D5_4ac3_BE4C_BF29F2E7AAE7

#include <vector>

#include "Node.h"
#include "Group.h"

class CNgroup:public CGroup{
public:
	CNgroup();
	virtual ~CNgroup();

private:
	vector<CNode*> mvNode;

public:
	void addNode(CNode* pNode);

	size_t getNumOfNode();
	CNode* getNode(size_t index);

};
#endif //include_guard


