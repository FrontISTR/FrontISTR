/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Quad2.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BA1ED7F1_AE6F_4c8a_90B2_3C7BFC9377D6
#define BA1ED7F1_AE6F_4c8a_90B2_3C7BFC9377D6

#include "Element.h"

class CQuad2 : public CElement
{
public:
	CQuad2(void);
	virtual ~CQuad2(void);

	virtual size_t getNodeID_Fistr2MW3(size_t index);

	virtual size_t getNumOfFace(){return 1;}
	virtual size_t getNumOfEdge(){ return 4;}
	virtual const char* getCType(){ return "Quad2";}
	virtual size_t getNType(){ return FistrElementType::Quad2;}

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










