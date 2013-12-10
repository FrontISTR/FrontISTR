/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   IFaceQuad2.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef AB7A00BD_1906_4d81_BEF7_4C07C5044F16
#define AB7A00BD_1906_4d81_BEF7_4C07C5044F16

#include "Quad2.h"

class CIFaceQuad2 : public CQuad2
{
public:
	CIFaceQuad2(void);
	virtual ~CIFaceQuad2(void);

	virtual const char* getCType(){ return "IFaceQuad2";}
	virtual size_t getNType(){ return FistrElementType::IFaceQuad2;}

	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard




















