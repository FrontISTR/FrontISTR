/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Triangle.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BFA325A7_5F85_4143_BD6F_D9A1D92B8188
#define BFA325A7_5F85_4143_BD6F_D9A1D92B8188

#include "Element.h"

class CTriangle : public CElement
{
public:
	CTriangle(void);
	virtual ~CTriangle(void);

	virtual size_t getNodeID_Fistr2MW3(size_t index);

	virtual size_t getNumOfFace(){return 1;}
	virtual size_t getNumOfEdge(){ return 3;}
	virtual const char* getCType(){ return "Triangle";}
	virtual size_t getNType(){ return FistrElementType::Triangle;}

	// Face
	virtual vector<CNode*> getFistrFaceNode(size_t nface);
	virtual size_t getMW3FaceNum(size_t fistr_nface);

	virtual size_t getFaceNodeNum(size_t fistr_nface);
	virtual string getFaceType(size_t fistr_nface);

	// Edge
	virtual vector<CNode*> getFistrEdgeNode(size_t nedge);
	virtual string getEdgeType(size_t nedge);
	virtual size_t getMW3EdgeNum(size_t fistr_nedge);
};
#endif // include guard













