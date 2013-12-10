/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   IFaceQuad.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef CE010_E927_4d38_8B27_AC76F7757CE1
#define CE010_E927_4d38_8B27_AC76F7757CE1

#include "Quad.h"

class CIFaceQuad : public CQuad
{
public:
	CIFaceQuad(void);
	virtual ~CIFaceQuad(void);

	virtual const char* getCType(){ return "IFaceQuad";}
	virtual size_t getNType(){ return FistrElementType::IFaceQuad;}

	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard














