/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommPyramid.h
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
#ifndef _COMMPYRAMID_H
#define	_COMMPYRAMID_H
class CCommPyramid:public CCommElement{
public:
    CCommPyramid();
    virtual ~CCommPyramid();
public:
    virtual bool isTypeCoincidence();
    virtual uint getShapeType(){ return ElementType::Pyramid;}
    virtual uint getBaseShapeType(){ return BaseElementType::Solid;}
};
#endif	/* _COMMPYRAMID_H */
}
