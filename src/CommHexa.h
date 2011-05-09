/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommHexa.h
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
#ifndef _COMMHEXA_H
#define	_COMMHEXA_H
class CCommHexa:public CCommElement{
public:
    CCommHexa();
    virtual ~CCommHexa();
public:
    virtual bool isTypeCoincidence();
    virtual uint getShapeType(){ return ElementType::Hexa;}
    virtual uint getBaseShapeType(){ return BaseElementType::Solid;}
};
#endif	/* _COMMHEXA_H */
}
