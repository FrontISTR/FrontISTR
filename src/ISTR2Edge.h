/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ISTR2Edge.h
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
#ifndef _ISTR2EDGE_H
#define	_ISTR2EDGE_H
class CISTR2Edge{
private:
    CISTR2Edge();
public:
    static CISTR2Edge* Instance(){
        static CISTR2Edge moISTR2Edge;
        return &moISTR2Edge;
    }
    ~CISTR2Edge();
private:
    Utility::CLogger *mpLogger;
    uint mvHexa[20];
    uint mvTetra[10];
    uint mvPrism[15];
    uint mvQuad[8];
    uint mvTriangle[6];
public:
    uint& HexaEdgeNum(const uint& ishape);
    uint& TetraEdgeNum(const uint& ishape);
    uint& PrismEdgeNum(const uint& ishape);
    uint& QuadEdgeNum(const uint& ishape);
    uint& TriangleEdgeNum(const uint& ishape);
};
#endif	/* _ISTR2EDGE_H */
}
