/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/CommTetra.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
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
#ifndef _COMMTETRA_H
#define	_COMMTETRA_H
class CCommTetra:public CCommElement{
public:
    CCommTetra();
    virtual ~CCommTetra();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType(){ return ElementType::Tetra;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Solid;}
};
#endif	/* _COMMTETRA_H */
}
