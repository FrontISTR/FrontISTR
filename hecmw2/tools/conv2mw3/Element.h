/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Element.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef A35E3DA_F1AE_40c8_AEE9_A616C7304E38
#define A35E3DA_F1AE_40c8_AEE9_A616C7304E38

#include <string>
#include <vector>
using namespace std;

#include <boost/lexical_cast.hpp>
using namespace boost;

#include "Node.h"
#include "FrontISTR_Type.h"

#include "Message.h"

class CElement{
public:
	CElement();
	virtual ~CElement();

protected:
	size_t mID;
	vector<CNode*> mvNode;

public:
	void setID(size_t id);
	size_t getID(){ return mID;}

	void setNode(size_t index, CNode* pNode);
	
	size_t getNumOfNode(){ return mvNode.size();}
	size_t getNodeID(size_t index){ return mvNode[index]->getID();}
	virtual size_t getNodeID_Fistr2MW3(size_t index)=0;

	virtual size_t getNumOfFace()=0;
	virtual size_t getNumOfEdge()=0;
	virtual const char* getCType()=0;
	virtual size_t getNType()=0;

	// Face
	virtual vector<CNode*> getFistrFaceNode(size_t nface)=0;
	virtual size_t getMW3FaceNum(size_t fistr_nface)=0;

	virtual size_t getFaceNodeNum(size_t fistr_nface)=0;
	virtual string getFaceType(size_t fistr_nface)=0;

	// Edge
	virtual vector<CNode*> getFistrEdgeNode(size_t nedge)=0;
	virtual string getEdgeType(size_t nedge)=0;
	virtual size_t getMW3EdgeNum(size_t fistr_nedge)=0;
};

#endif // include guard








