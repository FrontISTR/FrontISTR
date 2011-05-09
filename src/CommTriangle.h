/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommTriangle.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommElement.h"
namespace pmw{
#ifndef _COMMTRIANGLE_H
#define	_COMMTRIANGLE_H
class CCommTriangle:public CCommElement{
public:
    CCommTriangle();
    virtual ~CCommTriangle();
public:
    virtual bool isTypeCoincidence();
    virtual uint getShapeType(){ return ElementType::Triangle;}
    virtual uint getBaseShapeType(){ return BaseElementType::Shell;}
};
#endif	/* _COMMTRIANGLE_H */
}
