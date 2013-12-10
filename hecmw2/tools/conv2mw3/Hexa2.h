/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Hexa2.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef AA2D_D59C_457e_B0FE_4E60DEE9152F
#define AA2D_D59C_457e_B0FE_4E60DEE9152F

#include "Element.h"

class CHexa2 : public CElement
{
public:
	CHexa2(void);
	virtual ~CHexa2(void);

	virtual size_t getNodeID_Fistr2MW3(size_t index);

	virtual size_t getNumOfFace(){return 6;}
	virtual size_t getNumOfEdge(){ return 12;}
	virtual const char* getCType(){ return "Hexa2";}
	virtual size_t getNType(){ return FistrElementType::Hexa2;}

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















