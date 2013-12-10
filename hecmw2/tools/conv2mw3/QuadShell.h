/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   QuadShell.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef B1F33C82_342A_4fa4_B77E_11A128673394
#define B1F33C82_342A_4fa4_B77E_11A128673394

#include "Quad.h"

class CQuadShell : public CQuad
{
public:
	CQuadShell(void);
	virtual ~CQuadShell(void);

	virtual const char* getCType(){ return "QuadShell";}
	virtual size_t getNType(){ return FistrElementType::QuadShell;}

	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard









