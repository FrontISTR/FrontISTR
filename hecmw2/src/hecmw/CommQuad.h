/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommQuad.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommElement.h"
namespace pmw
{
#ifndef _COMMQUAD_H
#define	_COMMQUAD_H
class CCommQuad:public CCommElement
{
public:
    CCommQuad();
    virtual ~CCommQuad();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType() {
        return ElementType::Quad;
    }
    virtual uiint getBaseShapeType() {
        return BaseElementType::Shell;
    }
};
#endif	/* _COMMQUAD_H */
}
