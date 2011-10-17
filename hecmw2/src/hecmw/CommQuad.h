/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/CommQuad.h
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
namespace pmw{
#ifndef _COMMQUAD_H
#define	_COMMQUAD_H
class CCommQuad:public CCommElement{
public:
    CCommQuad();
    virtual ~CCommQuad();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType(){ return ElementType::Quad;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Shell;}
};
#endif	/* _COMMQUAD_H */
}
