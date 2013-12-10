/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   QuadShell2.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef E2DF39_88A1_4430_8EC7_DDC1AE189BC7
#define E2DF39_88A1_4430_8EC7_DDC1AE189BC7

#include "Quad2.h"

class CQuadShell2 : public CQuad2
{
public:
	CQuadShell2(void);
	virtual ~CQuadShell2(void);

	virtual const char* getCType(){ return "QuadShell2";}
	virtual size_t getNType(){ return FistrElementType::QuadShell2;}


	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard



