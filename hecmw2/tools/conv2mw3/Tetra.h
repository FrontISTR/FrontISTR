/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Tetra.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef CA5F3A_6493_4816_9E6E_CE6F25B7C6A9
#define CA5F3A_6493_4816_9E6E_CE6F25B7C6A9

#include "Element.h"

class CTetra : public CElement
{
public:
	CTetra(void);
	virtual ~CTetra(void);

	virtual size_t getNodeID_Fistr2MW3(size_t index);

	virtual size_t getNumOfFace(){return 4;}
	virtual size_t getNumOfEdge(){ return 6;}
	virtual const char* getCType(){ return "Tetra";}
	virtual size_t getNType(){ return FistrElementType::Tetra;}

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









