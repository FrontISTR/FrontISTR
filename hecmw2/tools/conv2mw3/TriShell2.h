/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   TriShell2.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BA9834_7CBF_4bbe_8C0A_F217AFBE14A4
#define BA9834_7CBF_4bbe_8C0A_F217AFBE14A4

#include "Triangle2.h"

class CTriShell2 : public CTriangle2
{
public:
	CTriShell2(void);
	virtual ~CTriShell2(void);

	virtual const char* getCType(){ return "TriShell2";}
	virtual size_t getNType(){ return FistrElementType::TriShell2;}


	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard









