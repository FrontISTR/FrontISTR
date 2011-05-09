
#include "Logger.h"

//
//  ISTR2Edge.cpp
//
//
//
//
//              2010.02.18
//              k.Takeda
#include "ISTR2Edge.h"
using namespace pmw;

CISTR2Edge::CISTR2Edge()
{
    uint i;
    // FrontISTR -> MW3 辺番号 (頂点Overの番号について辺番号)
    // --
    // Hexa 20 Node
    for(i=0; i< 20; i++){
        if(i< 8) mvHexa[i]= i;
        if(i>=8) mvHexa[i]= i-8;
    };
    // Tetra 10 Node
    for(i=0; i< 4; i++){
        mvTetra[i]= i;
    };
    mvTetra[4]= 1;
    mvTetra[5]= 2;
    mvTetra[6]= 0;
    mvTetra[7]= 3;
    mvTetra[8]= 4;
    mvTetra[9]= 5;
    // Prism 15 Node
    for(i=0; i< 6; i++){
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

    // Quad 8 Node
    for(i=0; i< 8; i++){
        if(i< 4) mvQuad[i]= i;
        if(i> 3) mvQuad[i]= i-4;
    };
    // Triangle 6 Node
    for(i=0; i< 6; i++){
        if(i< 3) mvTriangle[i]= i;
        if(i> 2) mvTriangle[i]= i-3;
    }

    mpLogger= Utility::CLogger::Instance();
}
CISTR2Edge::~CISTR2Edge()
{
    ;
}

// 形状関数 番号 -> MW3の辺番号に変換
// --
uint& CISTR2Edge::HexaEdgeNum(const uint& ishape)
{
    if(ishape < 8) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 7");
    if(ishape > 19){
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }
    
    return mvHexa[ishape];
}
uint& CISTR2Edge::TetraEdgeNum(const uint& ishape)
{
    if(ishape < 4) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 3");
    if(ishape > 9){
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }

    return mvTetra[ishape];
}
uint& CISTR2Edge::PrismEdgeNum(const uint& ishape)
{
    if(ishape < 6) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 5");
    if(ishape > 14){
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }

    return mvPrism[ishape];
}
uint& CISTR2Edge::QuadEdgeNum(const uint& ishape)
{
    if(ishape < 4) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 3");
    if(ishape > 7){
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }

    return mvQuad[ishape];
}
uint& CISTR2Edge::TriangleEdgeNum(const uint& ishape)
{
    if(ishape < 3) mpLogger->Info(Utility::LoggerMode::Warn,"edge is shape num > 2");
    if(ishape > 5){
        mpLogger->Info(Utility::LoggerMode::Error, "node size over");
        return mpLogger->getUDummyValue();
    }

    return mvTriangle[ishape];
}

















