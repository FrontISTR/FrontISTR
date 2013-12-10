/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   TriShell.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BCBD634_D05D_4e82_98EE_56EFC10A5AED
#define BCBD634_D05D_4e82_98EE_56EFC10A5AED

#include "Triangle.h"

class CTriShell : public CTriangle
{
public:
	CTriShell(void);
	virtual ~CTriShell(void);

	virtual const char* getCType(){ return "TriShell";}
	virtual size_t getNType(){ return FistrElementType::TriShell;}

	virtual string getFaceType(size_t fistr_nface);
};
#endif // include guard













