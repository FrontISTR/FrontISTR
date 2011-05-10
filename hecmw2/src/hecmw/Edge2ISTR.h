/* 
 * File:   Edge2ISTR.h
 * Author: ktakeda
 *
 * MW3 辺番号 -> FrontISTR 2次節点番号
 *
 * Created on 2010/02/18, 15:17
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

    // MW3 辺番号 -> FrontISTR 2次節点番号
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



