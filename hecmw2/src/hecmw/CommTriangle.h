/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommTriangle.h
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
#ifndef _COMMTRIANGLE_H
#define	_COMMTRIANGLE_H
class CCommTriangle:public CCommElement
{
public:
    CCommTriangle();
    virtual ~CCommTriangle();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType() {
        return ElementType::Triangle;
    }
    virtual uiint getBaseShapeType() {
        return BaseElementType::Shell;
    }
};
#endif	/* _COMMTRIANGLE_H */
}
