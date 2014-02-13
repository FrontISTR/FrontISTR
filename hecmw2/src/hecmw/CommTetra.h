/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommTetra.h
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
#include "CommHexa.h"
namespace pmw
{
#ifndef _COMMTETRA_H
#define	_COMMTETRA_H
class CCommTetra:public CCommElement
{
public:
    CCommTetra();
    virtual ~CCommTetra();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType() {
        return ElementType::Tetra;
    }
    virtual uiint getBaseShapeType() {
        return BaseElementType::Solid;
    }
};
#endif	/* _COMMTETRA_H */
}
