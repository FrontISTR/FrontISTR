//
//  ProgElementTree.cpp
//
//
//
//              2009.09.03
//              2009.09.03
//              k.Takeda
#include "ProgElementTree.h"
using namespace pmw;

// construct & destruct
//
CProgElementTree::CProgElementTree()
{
    uint ivert,iedge,iface,invalid;
    invalid=9;//progElem(Hexa)の”最大頂点数+1” をInvalid Value とする

    mInvalidNum= invalid;//prolongationの頂点が無い配列に使用.

    // == Vert ==

    // 子のAddress(親の頂点番号) => 子の頂点番号
    // --
    // Hexa
    // --
    mHexaVertChildVert[0]=0;
    mHexaVertChildVert[1]=1;
    mHexaVertChildVert[2]=2;
    mHexaVertChildVert[3]=3;
    mHexaVertChildVert[4]=4;
    mHexaVertChildVert[5]=5;
    mHexaVertChildVert[6]=6;
    mHexaVertChildVert[7]=7;
    
    // --
    // Tetra
    // --
    mTetraVertChildVert[0]=1;
    mTetraVertChildVert[1]=2;
    mTetraVertChildVert[2]=0;
    mTetraVertChildVert[3]=5;
    
    // --
    // Prism
    // --
    mPrismVertChildVert[0]=1;
    mPrismVertChildVert[1]=2;
    mPrismVertChildVert[2]=0;
    mPrismVertChildVert[3]=5;
    mPrismVertChildVert[4]=6;
    mPrismVertChildVert[5]=4;
    
    // --
    // Pyramid
    // --
    mPyramidVertChildVert[0]=0;
    mPyramidVertChildVert[1]=1;
    mPyramidVertChildVert[2]=2;
    mPyramidVertChildVert[3]=3;
    mPyramidVertChildVert[4]=1;
    
    // --
    // Quad
    // --
    mQuadVertChildVert[0]=0;
    mQuadVertChildVert[1]=1;
    mQuadVertChildVert[2]=1;
    mQuadVertChildVert[3]=1;
    
    // --
    // Triangle
    // --
    mTriangleVertChildVert[0]=1;
    mTriangleVertChildVert[1]=1;
    mTriangleVertChildVert[2]=1;
    
    // --
    // Beam
    // --
    mBeamVertChildVert[0]=0;
    mBeamVertChildVert[1]=1;
    

    // == Edge ==
    
    // --
    // Hexa - Edge
    // --

    // 子(progElem)のAddress(親の頂点番号)-親の辺番号 => 子(progElem)の頂点番号
    // --
    for(ivert=0; ivert< 8; ivert++){
        for(iedge=0; iedge< 12; iedge++){
            mHexaEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    // 辺番号に対応する,子の頂点番号
    // --
    // 例:親の辺番号0,3,8に相当する,子の頂点番号 : 親の頂点0に居る子(progElem)
    // --
    // 子(progElem)は,8個のHexa
    mHexaEdgeChildVert[0][0]=1; mHexaEdgeChildVert[0][3]=3; mHexaEdgeChildVert[0][8]=4;//子の住所:親の頂点0
    mHexaEdgeChildVert[1][0]=0; mHexaEdgeChildVert[1][1]=2; mHexaEdgeChildVert[1][9]=5;
    mHexaEdgeChildVert[2][1]=1; mHexaEdgeChildVert[2][2]=3; mHexaEdgeChildVert[2][10]=6;
    mHexaEdgeChildVert[3][3]=0; mHexaEdgeChildVert[3][2]=2; mHexaEdgeChildVert[3][11]=7;
    mHexaEdgeChildVert[4][8]=0; mHexaEdgeChildVert[4][4]=5; mHexaEdgeChildVert[4][7]=7;
    mHexaEdgeChildVert[5][9]=1; mHexaEdgeChildVert[5][4]=4; mHexaEdgeChildVert[5][5]=6;
    mHexaEdgeChildVert[6][10]=2; mHexaEdgeChildVert[6][5]=5; mHexaEdgeChildVert[6][6]=7;
    mHexaEdgeChildVert[7][11]=3; mHexaEdgeChildVert[7][7]=4; mHexaEdgeChildVert[7][6]=6;
    

    // --
    // Tetra - Edge
    // --
    for(ivert=0; ivert< 4; ivert++){
        for(iedge=0; iedge< 6; iedge++){
            mTetraEdgeChildVert[ivert][iedge]= invalid;
        }
    }
    //子(progElem)は,4個のHexa
    mTetraEdgeChildVert[0][2]=0; mTetraEdgeChildVert[0][0]=2; mTetraEdgeChildVert[0][3]=5;
    mTetraEdgeChildVert[1][0]=1; mTetraEdgeChildVert[1][1]=3; mTetraEdgeChildVert[1][4]=6;
    mTetraEdgeChildVert[2][2]=1; mTetraEdgeChildVert[2][1]=3; mTetraEdgeChildVert[2][5]=4;
    mTetraEdgeChildVert[3][3]=1; mTetraEdgeChildVert[3][5]=4; mTetraEdgeChildVert[3][4]=6;
    

    // --
    // Prism - Edge
    // --
    for(ivert=0; ivert< 6; ivert++){
        for(iedge=0; iedge< 9; iedge++){
            mPrismEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    //子(progElem)は,6個のHexa
    mPrismEdgeChildVert[0][1]=0; mPrismEdgeChildVert[0][0]=2; mPrismEdgeChildVert[0][3]=5;
    mPrismEdgeChildVert[1][0]=1; mPrismEdgeChildVert[1][2]=3; mPrismEdgeChildVert[1][4]=6;
    mPrismEdgeChildVert[2][1]=1; mPrismEdgeChildVert[2][2]=3; mPrismEdgeChildVert[2][5]=4;
    mPrismEdgeChildVert[3][3]=1; mPrismEdgeChildVert[3][8]=4; mPrismEdgeChildVert[3][6]=6;
    mPrismEdgeChildVert[4][4]=2; mPrismEdgeChildVert[4][6]=5; mPrismEdgeChildVert[4][7]=7;
    mPrismEdgeChildVert[5][5]=0; mPrismEdgeChildVert[5][8]=5; mPrismEdgeChildVert[5][7]=7;

    
    // --
    // Pyramid - Edge
    // --
    for(ivert=0; ivert< 5; ivert++){
        for(iedge=0; iedge< 8; iedge++){
            mPyramidEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    //子(progElem)は,4個のHexa && 4個のPyramid
    mPyramidEdgeChildVert[0][0]=1; mPyramidEdgeChildVert[0][3]=3; mPyramidEdgeChildVert[0][7]=4;
    mPyramidEdgeChildVert[1][0]=0; mPyramidEdgeChildVert[1][1]=2; mPyramidEdgeChildVert[1][4]=5;
    mPyramidEdgeChildVert[2][1]=1; mPyramidEdgeChildVert[2][2]=3; mPyramidEdgeChildVert[2][5]=6;
    mPyramidEdgeChildVert[3][3]=0; mPyramidEdgeChildVert[3][2]=2; mPyramidEdgeChildVert[3][6]=7;
    mPyramidEdgeChildVert[4][4]=0; mPyramidEdgeChildVert[4][5]=2;
    mPyramidEdgeChildVert[5][5]=0; mPyramidEdgeChildVert[5][6]=2;
    mPyramidEdgeChildVert[6][6]=0; mPyramidEdgeChildVert[6][7]=2;
    mPyramidEdgeChildVert[7][7]=0; mPyramidEdgeChildVert[7][4]=2;

    
    // --
    // Quad - Edge
    // --
    for(ivert=0; ivert< 4; ivert++){
        for(iedge=0; iedge< 4; iedge++){
            mQuadEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    mQuadEdgeChildVert[0][0]=1; mQuadEdgeChildVert[0][3]=3;
    mQuadEdgeChildVert[1][0]=0; mQuadEdgeChildVert[1][1]=2;
    mQuadEdgeChildVert[2][1]=0; mQuadEdgeChildVert[2][2]=2;
    mQuadEdgeChildVert[3][2]=0; mQuadEdgeChildVert[3][3]=2;
    
    
    // --
    // Triangle - Edge
    // --
    for(ivert=0; ivert< 3; ivert++){
        for(iedge=0; iedge< 3; iedge++){
            mTriangleEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    mTriangleEdgeChildVert[0][2]=0; mTriangleEdgeChildVert[0][0]=2;
    mTriangleEdgeChildVert[1][0]=0; mTriangleEdgeChildVert[1][1]=2;
    mTriangleEdgeChildVert[2][1]=0; mTriangleEdgeChildVert[2][2]=2;
    
    
    // --
    // Beam - Edge
    // --
    for(ivert=0; ivert< 2; ivert++){
        for(iedge=0; iedge< 1; iedge++){
            mBeamEdgeChildVert[ivert][iedge]=invalid;
        }
    }
    mBeamEdgeChildVert[0][0]=1;
    mBeamEdgeChildVert[1][0]=0;


    // == Face ==
    
    // --
    // Hexa - Face
    // --

    // 子(progElem)のAddress(親の頂点番号)-親の面番号 => 子(progElem)の頂点番号
    // --
    for(ivert=0; ivert< 8; ivert++){
        for(iface=0; iface< 6; iface++){
            mHexaFaceChildVert[ivert][iface]=invalid;
        }
    }
    // 面番号に対応する,子の頂点番号
    // --
    // 例:親の面番号0,4,3に相当する,子の頂点番号 : 親の頂点0に居る子(progElem)
    // {データ入力が,頂点順になっていないのは,Factoryで適当に生成順を決めてしまったから.}
    // --
    mHexaFaceChildVert[0][0]=2; mHexaFaceChildVert[0][4]=5; mHexaFaceChildVert[0][3]=7;
    mHexaFaceChildVert[1][0]=3; mHexaFaceChildVert[1][4]=4; mHexaFaceChildVert[1][2]=6;
    mHexaFaceChildVert[2][0]=0; mHexaFaceChildVert[2][2]=5; mHexaFaceChildVert[2][5]=7;
    mHexaFaceChildVert[3][0]=1; mHexaFaceChildVert[3][3]=4; mHexaFaceChildVert[3][5]=6;
    mHexaFaceChildVert[4][4]=1; mHexaFaceChildVert[4][3]=3; mHexaFaceChildVert[4][1]=6;
    mHexaFaceChildVert[5][4]=0; mHexaFaceChildVert[5][2]=2; mHexaFaceChildVert[5][1]=7;
    mHexaFaceChildVert[6][2]=1; mHexaFaceChildVert[6][5]=3; mHexaFaceChildVert[6][1]=4;
    mHexaFaceChildVert[7][3]=0; mHexaFaceChildVert[7][5]=2; mHexaFaceChildVert[7][1]=5;
    

    // --
    // Tetra - Face
    // --
    for(ivert=0; ivert< 4; ivert++){
        for(iface=0; iface< 4; iface++){
            mTetraFaceChildVert[ivert][iface]=invalid;
        }
    }
    mTetraFaceChildVert[0][0]=3; mTetraFaceChildVert[0][3]=4; mTetraFaceChildVert[0][1]=6;
    mTetraFaceChildVert[1][0]=0; mTetraFaceChildVert[1][1]=5; mTetraFaceChildVert[1][2]=7;
    mTetraFaceChildVert[2][0]=2; mTetraFaceChildVert[2][3]=5; mTetraFaceChildVert[2][2]=7;
    mTetraFaceChildVert[3][3]=0; mTetraFaceChildVert[3][1]=2; mTetraFaceChildVert[3][2]=7;

    // --
    // Prism - Face
    // --
    for(ivert=0; ivert< 6; ivert++){
        for(iface=0; iface< 5; iface++){
            mPrismFaceChildVert[ivert][iface]=invalid;
        }
    }
    mPrismFaceChildVert[0][0]=3; mPrismFaceChildVert[0][4]=4; mPrismFaceChildVert[0][2]=6;
    mPrismFaceChildVert[1][0]=0; mPrismFaceChildVert[1][2]=5; mPrismFaceChildVert[1][3]=7;
    mPrismFaceChildVert[2][0]=2; mPrismFaceChildVert[2][4]=5; mPrismFaceChildVert[2][3]=7;
    mPrismFaceChildVert[3][4]=0; mPrismFaceChildVert[3][2]=2; mPrismFaceChildVert[3][1]=7;
    mPrismFaceChildVert[4][2]=1; mPrismFaceChildVert[4][3]=3; mPrismFaceChildVert[4][1]=4;
    mPrismFaceChildVert[5][4]=1; mPrismFaceChildVert[5][3]=3; mPrismFaceChildVert[5][1]=6;

    // --
    // Pyramid - Face
    // --
    for(ivert=0; ivert< 5; ivert++){
        for(iface=0; iface< 5; iface++){
            mPyramidFaceChildVert[ivert][iface]=invalid;
        }
    }
    mPyramidFaceChildVert[0][0]=2; mPyramidFaceChildVert[0][4]=5; mPyramidFaceChildVert[0][3]=7;
    mPyramidFaceChildVert[1][0]=3; mPyramidFaceChildVert[1][4]=4; mPyramidFaceChildVert[1][1]=6;
    mPyramidFaceChildVert[2][0]=0; mPyramidFaceChildVert[2][1]=5; mPyramidFaceChildVert[2][2]=7;
    mPyramidFaceChildVert[3][0]=1; mPyramidFaceChildVert[3][3]=4; mPyramidFaceChildVert[3][2]=6;
    mPyramidFaceChildVert[4][1]=3;
    mPyramidFaceChildVert[5][2]=3;
    mPyramidFaceChildVert[6][3]=3;
    mPyramidFaceChildVert[7][4]=3;
   
    // --
    // Quad - Face
    // --
    for(ivert=0; ivert< 4; ivert++){
        for(iface=0; iface< 1; iface++){
            mQuadFaceChildVert[ivert][iface]=invalid;
        }
    }
    mQuadFaceChildVert[0][0]=2;
    mQuadFaceChildVert[1][0]=3;
    mQuadFaceChildVert[2][0]=3;
    mQuadFaceChildVert[3][0]=3;
    
    // --
    // Triangle - Face
    // --
    for(ivert=0; ivert< 3; ivert++){
        for(iface=0; iface< 1; iface++){
            mTriangleFaceChildVert[ivert][iface]=invalid;
        }
    }
    mTriangleFaceChildVert[0][0]=3;
    mTriangleFaceChildVert[1][0]=3;
    mTriangleFaceChildVert[2][0]=3;
    

    // == Vol Point ==
    
    // --
    // Hexa Vol-Point
    // --
    mHexaVolChildVert[0]=6;
    mHexaVolChildVert[1]=7;
    mHexaVolChildVert[2]=4;
    mHexaVolChildVert[3]=5;
    mHexaVolChildVert[4]=2;
    mHexaVolChildVert[5]=3;
    mHexaVolChildVert[6]=0;
    mHexaVolChildVert[7]=1;

    // --
    // Tetra Vol-Point
    // --
    mTetraVolChildVert[0]=7;
    mTetraVolChildVert[1]=4;
    mTetraVolChildVert[2]=6;
    mTetraVolChildVert[3]=3;

    // --
    // Prism Vol-Point
    // --
    mPrismVolChildVert[0]=7;
    mPrismVolChildVert[1]=4;
    mPrismVolChildVert[2]=6;
    mPrismVolChildVert[3]=3;
    mPrismVolChildVert[4]=0;
    mPrismVolChildVert[5]=2;

    // --
    // Pyramid Vol-Point
    // --
    mPyramidVolChildVert[0]=6;
    mPyramidVolChildVert[1]=7;
    mPyramidVolChildVert[2]=4;
    mPyramidVolChildVert[3]=5;
    mPyramidVolChildVert[4]=4;
    mPyramidVolChildVert[5]=4;
    mPyramidVolChildVert[6]=4;
    mPyramidVolChildVert[7]=4;
    
    // --
    // Quad Vol-Point
    // --
    mQuadVolChildVert[0]=2;
    mQuadVolChildVert[1]=3;
    mQuadVolChildVert[2]=3;
    mQuadVolChildVert[3]=3;

    // --
    // Triangle Vol-Point
    // --
    mTriangleVolChildVert[0]=3;
    mTriangleVolChildVert[1]=3;
    mTriangleVolChildVert[2]=3;

    // --
    // Beam Vol-Point
    // --
    mBeamVolChildVert[0]=1;
    mBeamVolChildVert[1]=0;
    
}

CProgElementTree::~CProgElementTree()
{
    ;
}

// 親の頂点番号=> 子の頂点番号
uint& CProgElementTree::getVertProgVert(const uint& ivert, const uint& elemType)
{
    switch(elemType){
        case(ElementType::Hexa):
            mProgVert= mHexaVertChildVert[ivert];
            break;
        case(ElementType::Tetra):
            mProgVert= mTetraVertChildVert[ivert];
            break;
        case(ElementType::Prism):
            mProgVert= mPrismVertChildVert[ivert];
            break;
        case(ElementType::Pyramid):
            mProgVert= mPyramidVertChildVert[ivert];
            break;
        case(ElementType::Quad):
            mProgVert= mQuadVertChildVert[ivert];
            break;
        case(ElementType::Triangle):
            mProgVert= mTriangleVertChildVert[ivert];
            break;
        case(ElementType::Beam):
            mProgVert= mBeamVertChildVert[ivert];
            break;
        default:
            break;
    }
    return mProgVert;
}


// 親の辺番号と子のアドレス(親の頂点番号) => 子供の頂点番号
uint& CProgElementTree::getEdgeProgVert(const uint& iedge, const uint& child_address, const uint& elemType)
{
    switch(elemType){
        case(ElementType::Hexa):
            mProgVert= mHexaEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Tetra):
            mProgVert= mTetraEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Prism):
            mProgVert= mPrismEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Pyramid):
            mProgVert= mPyramidEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Quad):
            mProgVert= mQuadEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Triangle):
            mProgVert= mTriangleEdgeChildVert[child_address][iedge];
            break;
        case(ElementType::Beam):
            mProgVert= mBeamEdgeChildVert[child_address][iedge];
            break;
        default:
            break;
    }

    return mProgVert;
}

uint& CProgElementTree::getFaceProgVert(const uint& iface, const uint& child_address, const uint& elemType)
{
    Utility::CLogger *pLogger= Utility::CLogger::Instance();
    
    switch(elemType){
        case(ElementType::Hexa):
            mProgVert= mHexaFaceChildVert[child_address][iface];
            break;
        case(ElementType::Tetra):
            mProgVert= mTetraFaceChildVert[child_address][iface];
            break;
        case(ElementType::Prism):
            mProgVert= mPrismFaceChildVert[child_address][iface];
            break;
        case(ElementType::Pyramid):
            mProgVert= mPyramidFaceChildVert[child_address][iface];
            break;
        case(ElementType::Quad):
            mProgVert= mQuadFaceChildVert[child_address][iface];
            break;
        case(ElementType::Triangle):
            mProgVert= mTriangleFaceChildVert[child_address][iface];
            break;
        case(ElementType::Beam):
            //Error
            pLogger->Info(Utility::LoggerMode::Error, "at CProgElementTree::getFaceProgVert");
            break;
        default:
            break;
    }

    return mProgVert;
}

uint& CProgElementTree::getVolProgVert(const uint& child_address, const uint& elemType)
{
    switch(elemType){
        case(ElementType::Hexa):
            mProgVert = mHexaVolChildVert[child_address];
            break;
        case(ElementType::Tetra):
            mProgVert = mTetraVolChildVert[child_address];
            break;
        case(ElementType::Prism):
            mProgVert = mPrismVolChildVert[child_address];
            break;
        case(ElementType::Pyramid):
            mProgVert = mPyramidVolChildVert[child_address];
            break;
        case(ElementType::Quad):
            mProgVert = mQuadVolChildVert[child_address];
            break;
        case(ElementType::Triangle):
            mProgVert = mTriangleVolChildVert[child_address];
            break;
        case(ElementType::Beam):
            mProgVert = mBeamVolChildVert[child_address];
            break;
        default:
            break;
    }

    return mProgVert;
}


























