/* 
 * File:   ISTR2Edge.h
 * Author: ktakeda
 *
 * 節点番号の変換
 *　FrontISTRの節点番号(形状関数のShape番号) => MW3の頂点番号,辺番号
 *
 * Created on 2010/02/18, 15:14
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

    // FrontISTR -> MW3 辺番号 (頂点Overの番号について辺番号)
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


