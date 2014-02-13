/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ISTR2Edge.h
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
#ifndef _ISTR2EDGE_H
#define	_ISTR2EDGE_H
class CISTR2Edge
{
private:
    CISTR2Edge();
public:
    static CISTR2Edge* Instance() {
        static CISTR2Edge moISTR2Edge;
        return &moISTR2Edge;
    }
    ~CISTR2Edge();
private:
    Utility::CLogger *mpLogger;
    uiint mvHexa[20];
    uiint mvTetra[10];
    uiint mvPrism[15];
    uiint mvQuad[8];
    uiint mvTriangle[6];
public:
    uiint& HexaEdgeNum(const uiint& ishape);
    uiint& TetraEdgeNum(const uiint& ishape);
    uiint& PrismEdgeNum(const uiint& ishape);
    uiint& QuadEdgeNum(const uiint& ishape);
    uiint& TriangleEdgeNum(const uiint& ishape);
};
#endif	/* _ISTR2EDGE_H */
}
