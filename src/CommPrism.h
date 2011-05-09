/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommPrism.h
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
#include "CommHexa.h"
namespace pmw{
#ifndef _COMMPRISM_H
#define	_COMMPRISM_H
class CCommPrism:public CCommElement{
public:
    CCommPrism();
    virtual ~CCommPrism();
public:
    virtual bool isTypeCoincidence();
    virtual uint getShapeType(){ return ElementType::Prism;}
    virtual uint getBaseShapeType(){ return BaseElementType::Solid;}
};
#endif	/* _COMMPRISM_H */
}
