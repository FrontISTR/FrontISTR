/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Edge2ISTR.h
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
#include "TypeDef.h"
#include "Logger.h"
namespace pmw
{
#ifndef _EDGE2ISTR_H
#define	_EDGE2ISTR_H
class CEdge2ISTR
{
private:
    CEdge2ISTR();
public:
    static CEdge2ISTR* Instance() {
        static CEdge2ISTR moEdge2ISTR;
        return &moEdge2ISTR;
    }
    ~CEdge2ISTR();
private:
    Utility::CLogger *mpLogger;
    uiint mvHexa[12];
    uiint mvTetra[6];
    uiint mvPrism[9];
    uiint mvQuad[4];
    uiint mvTriangle[3];
public:
    uiint& HexaShapeNum(const uiint& iedge);
    uiint& TetraShapeNum(const uiint& iedge);
    uiint& PrismShapeNum(const uiint& iedge);
    uiint& QuadShapeNum(const uiint& iedge);
    uiint& TriangleShapeNum(const uiint& iedge);
};
#endif	/* _EDGE2ISTR_H */
}
