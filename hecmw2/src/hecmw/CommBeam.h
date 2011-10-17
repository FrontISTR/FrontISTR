/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/CommBeam.h
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
#ifndef _COMMBEAM_H
#define	_COMMBEAM_H
class CCommBeam:public CCommElement{
public:
    CCommBeam();
    virtual ~CCommBeam();
public:
    virtual bool isTypeCoincidence();
    virtual uiint getShapeType(){ return ElementType::Beam;}
    virtual uiint getBaseShapeType(){ return BaseElementType::Beam;}
};
#endif	/* _COMMBEAM_H */
}
