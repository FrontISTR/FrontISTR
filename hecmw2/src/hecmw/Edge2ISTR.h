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



