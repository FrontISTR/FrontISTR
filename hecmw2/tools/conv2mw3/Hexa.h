/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Hexa.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef A7D37C9_7E5A_4340_8210_AE103DD9108F
#define A7D37C9_7E5A_4340_8210_AE103DD9108F

#include "Element.h"

class CHexa : public CElement
{
public:
	CHexa(void);
	virtual ~CHexa(void);

	virtual size_t getNodeID_Fistr2MW3(size_t index);

	virtual size_t getNumOfFace(){return 6;}
	virtual size_t getNumOfEdge(){ return 12;}
	virtual const char* getCType(){ return "Hexa";}
	virtual size_t getNType(){ return FistrElementType::Hexa;}

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












