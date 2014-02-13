/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommPrism.h
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
#ifndef _COMMPRISM_H
#define	_COMMPRISM_H
class CCommPrism:public CCommElement
{
public:
    CCommPrism();
    virtual ~CCommPrism();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType() {
        return ElementType::Prism;
    }
    virtual uiint getBaseShapeType() {
        return BaseElementType::Solid;
    }
};
#endif	/* _COMMPRISM_H */
}
