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


