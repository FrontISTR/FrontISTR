/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommHexa.h
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
#ifndef _COMMHEXA_H
#define	_COMMHEXA_H
class CCommHexa:public CCommElement
{
public:
    CCommHexa();
    virtual ~CCommHexa();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType() {
        return ElementType::Hexa;
    }
    virtual uiint getBaseShapeType() {
        return BaseElementType::Solid;
    }
};
#endif	/* _COMMHEXA_H */
}
