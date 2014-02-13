/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ISTR2Edge.cpp
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
#include "Logger.h"
#include "ISTR2Edge.h"
using namespace pmw;
CISTR2Edge::CISTR2Edge()
{
    uiint i;
    for(i=0; i< 20; i++) {
        if(i< 8) mvHexa[i]= i;
        if(i>=8) mvHexa[i]= i-8;
    };
    for(i=0; i< 4; i++) {
        mvTetra[i]= i;
    };
    mvTetra[4]= 1;
    mvTetra[5]= 2;
    mvTetra[6]= 0;
    mvTetra[7]= 3;
    mvTetra[8]= 4;
    mvTetra[9]= 5;
    for(i=0; i< 6; i++) {
        mvPrism[i]= i;
    };
    mvPrism[6]=  2;
    mvPrism[7]=  1;
    mvPrism[8]=  0;
    mvPrism[9]=  7;
    mvPrism[10]= 8;
    mvPrism[11]= 6;
    mvPrism[12]= 3;
    mvPrism[13]= 4;
    mvPrism[14]= 5;
    for(i=0; i< 8; i++) {
        if(i< 4) mvQuad[i]= i;
        if(i> 3) mvQuad[i]= i-4;
    };
    for(i=0; i< 6; i++) {
        if(i< 3) mvTriangle[i]= i;
        if(i> 2) mvTriangle[i]= i-3;
    }
    mpLogger= Utility::CLogger::Instance();
}
CISTR2Edge::~CISTR2Edge()
{
    ;
}
uiint& CISTR2Edge::HexaEdgeNum(const uiint& ishape)
{
    if(ishape < 8) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 7");
    if(ishape > 19) {
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    return mvHexa[ishape];
}
uiint& CISTR2Edge::TetraEdgeNum(const uiint& ishape)
{
    if(ishape < 4) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 3");
    if(ishape > 9) {
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    return mvTetra[ishape];
}
uiint& CISTR2Edge::PrismEdgeNum(const uiint& ishape)
{
    if(ishape < 6) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 5");
    if(ishape > 14) {
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    return mvPrism[ishape];
}
uiint& CISTR2Edge::QuadEdgeNum(const uiint& ishape)
{
    if(ishape < 4) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 3");
    if(ishape > 7) {
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    return mvQuad[ishape];
}
uiint& CISTR2Edge::TriangleEdgeNum(const uiint& ishape)
{
    if(ishape < 3) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 2");
    if(ishape > 5) {
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    return mvTriangle[ishape];
}
