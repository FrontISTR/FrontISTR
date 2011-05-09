/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Edge2ISTR.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "Logger.h"
namespace pmw{
#ifndef _EDGE2ISTR_H
#define	_EDGE2ISTR_H
class CEdge2ISTR{
private:
    CEdge2ISTR();
public:
    static CEdge2ISTR* Instance(){
        static CEdge2ISTR moEdge2ISTR;
        return &moEdge2ISTR;
    }
    ~CEdge2ISTR();
private:
    Utility::CLogger *mpLogger;
    uint mvHexa[12];
    uint mvTetra[6];
    uint mvPrism[9];
    uint mvQuad[4];
    uint mvTriangle[3];
public:
    uint& HexaShapeNum(const uint& iedge);
    uint& TetraShapeNum(const uint& iedge);
    uint& PrismShapeNum(const uint& iedge);
    uint& QuadShapeNum(const uint& iedge);
    uint& TriangleShapeNum(const uint& iedge);
};
#endif	/* _EDGE2ISTR_H */
}
