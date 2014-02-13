/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Edge2ISTR.cpp
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
#include "Edge2ISTR.h"
using namespace pmw;
CEdge2ISTR::CEdge2ISTR()
{
    mvHexa[0]=   8;
    mvHexa[1]=   9;
    mvHexa[2]=  10;
    mvHexa[3]=  11;
    mvHexa[4]=  12;
    mvHexa[5]=  13;
    mvHexa[6]=  14;
    mvHexa[7]=  15;
    mvHexa[8]=  16;
    mvHexa[9]=  17;
    mvHexa[10]= 18;
    mvHexa[11]= 19;
    mvTetra[0]= 6;
    mvTetra[1]= 4;
    mvTetra[2]= 5;
    mvTetra[3]= 7;
    mvTetra[4]= 8;
    mvTetra[5]= 9;
    mvPrism[0]=  8;
    mvPrism[1]=  7;
    mvPrism[2]=  6;
    mvPrism[3]= 12;
    mvPrism[4]= 13;
    mvPrism[5]= 14;
    mvPrism[6]= 11;
    mvPrism[7]=  9;
    mvPrism[8]= 10;
    mvQuad[0]= 4;
    mvQuad[1]= 5;
    mvQuad[2]= 6;
    mvQuad[3]= 7;
    mvTriangle[0]= 3;
    mvTriangle[1]= 4;
    mvTriangle[2]= 5;
    mpLogger= Utility::CLogger::Instance();
}
CEdge2ISTR::~CEdge2ISTR()
{
    ;
}
uiint& CEdge2ISTR::HexaShapeNum(const uiint& iedge)
{
    if(iedge > 11) {
        mpLogger->Info(Utility::LoggerMode::Error,"edge size over, edge num < 12");
        return mpLogger->getUDummyValue();
    }
    return mvHexa[iedge];
}
uiint& CEdge2ISTR::TetraShapeNum(const uiint& iedge)
{
    if(iedge > 5) {
        mpLogger->Info(Utility::LoggerMode::Error,"edge size over, edge num < 6");
        return mpLogger->getUDummyValue();
    }
    return mvTetra[iedge];
}
uiint& CEdge2ISTR::PrismShapeNum(const uiint& iedge)
{
    if(iedge > 8) {
        mpLogger->Info(Utility::LoggerMode::Error,"edge size over, edge num < 9");
        return mpLogger->getUDummyValue();
    }
    return mvPrism[iedge];
}
uiint& CEdge2ISTR::QuadShapeNum(const uiint& iedge)
{
    if(iedge > 3) {
        mpLogger->Info(Utility::LoggerMode::Error,"edge size over, edge num < 4");
        return mpLogger->getUDummyValue();
    }
    return mvQuad[iedge];
}
uiint& CEdge2ISTR::TriangleShapeNum(const uiint& iedge)
{
    if(iedge > 2) {
        mpLogger->Info(Utility::LoggerMode::Error,"edge size over, edge num < 3");
        return mpLogger->getUDummyValue();
    }
    return mvTriangle[iedge];
}
