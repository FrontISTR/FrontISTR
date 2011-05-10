//
//  ShapeHexa.cpp
//
//  Hexa要素の形状関数
//
//
//                  2010.01.21
//                  k.Takeda
#include <math.h>
#include <vector>

#include "ShapeHexa.h"
using namespace pmw;

CShapeHexa::CShapeHexa()
{
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    mvEdge_2_ISTR[0]=   8;
//    mvEdge_2_ISTR[1]=   9;
//    mvEdge_2_ISTR[2]=  10;
//    mvEdge_2_ISTR[3]=  11;
//    mvEdge_2_ISTR[4]=  12;
//    mvEdge_2_ISTR[5]=  13;
//    mvEdge_2_ISTR[6]=  14;
//    mvEdge_2_ISTR[7]=  15;
//    mvEdge_2_ISTR[8]=  16;
//    mvEdge_2_ISTR[9]=  17;
//    mvEdge_2_ISTR[10]= 18;
//    mvEdge_2_ISTR[11]= 19;
//
//    // FrontISTR -> MW3 辺番号
//    uint i;
//    for(i=0; i< 20; i++){
//        if(i< 8) mvISTR_2_Edge[i]= i;
//        if(i>=8) mvISTR_2_Edge[i]= i-8;
//    };
    
    uint numOfIntg, numOfShape, dof;
    // vector領域確保
    // --
    // 8節点要素 積分点1
    // --
    numOfIntg= 1; numOfShape= 8; dof= 3;
    ResizeShape(mvN81, numOfIntg, numOfShape);   
    ResizeDeriv(mvdNdr81, numOfIntg, numOfShape, dof);
    // 8節点要素 積分点8
    // --
    numOfIntg= 8; numOfShape= 8; dof= 3;
    ResizeShape(mvN82, numOfIntg, numOfShape);   
    ResizeDeriv(mvdNdr82, numOfIntg, numOfShape, dof);
    
    // 20節点要素 積分点1
    // --
    numOfIntg= 1; numOfShape= 20; dof= 3;
    ResizeShape(mvN201,  numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr201, numOfIntg, numOfShape, dof);
    // 20節点要素 積分点8(1次元では2)
    //  --
    numOfIntg= 8; numOfShape= 20; dof= 3;
    ResizeShape(mvN202,  numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr202, numOfIntg, numOfShape, dof);
    // 20節点要素 積分点27(1次元では3)
    // --
    numOfIntg= 27; numOfShape= 20; dof= 3;
    ResizeShape(mvN203, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr203, numOfIntg, numOfShape, dof);
    
    //ヤコビアン
    mpJacobi81= new CJacobian;
    mpJacobi82= new CJacobian;
    mpJacobi201= new CJacobian;
    mpJacobi202= new CJacobian;
    mpJacobi203= new CJacobian;

    mpJacobi81->setupRegion(1,8);
    mpJacobi82->setupRegion(8,8);
    mpJacobi201->setupRegion(1,20);
    mpJacobi202->setupRegion(8,20);
    mpJacobi203->setupRegion(27,20);
    //--
    //空間座標(x,y,z)の導関数配列確保
    //--
    // 8節点要素 積分点1
    numOfIntg= 1; numOfShape= 8; dof= 3;
    ResizeDeriv(mvdNdx81, numOfIntg, numOfShape, dof);
    // 8節点要素 積分点8
    numOfIntg= 8; numOfShape= 8; dof= 3;
    ResizeDeriv(mvdNdx82, numOfIntg, numOfShape, dof);
    // 20節点要素 積分点1
    numOfIntg= 1; numOfShape= 20; dof= 3;
    ResizeDeriv(mvdNdx201, numOfIntg, numOfShape, dof);
    // 20節点要素 積分点8
    numOfIntg= 8; numOfShape= 20; dof= 3;
    ResizeDeriv(mvdNdx202, numOfIntg, numOfShape, dof);
    // 20節点要素 積分点27
    numOfIntg= 27; numOfShape= 20; dof= 3;
    ResizeDeriv(mvdNdx203, numOfIntg, numOfShape, dof);

    //--
    // detJの領域確保
    //--
    numOfIntg= 1; numOfShape= 8;
    Resize_detJ(mv_detJ81, numOfIntg);
    numOfIntg= 8; numOfShape= 8;
    Resize_detJ(mv_detJ82, numOfIntg);
    numOfIntg= 1; numOfShape= 20;
    Resize_detJ(mv_detJ201, numOfIntg);
    numOfIntg= 8; numOfShape= 20;
    Resize_detJ(mv_detJ202, numOfIntg);
    numOfIntg= 27; numOfShape= 20;
    Resize_detJ(mv_detJ203, numOfIntg);

    
    //積分座標(自然座標,正規化座標):1次元
    mGzi1[0]= 0.0;

    mGzi2[0]= -0.5773502691;
    mGzi2[1]=  0.5773502691;
    
    mGzi3[0]= -0.77459666923;
    mGzi3[1]=  0.0000000000;
    mGzi3[2]=  0.77459666923;
    
    mGzi4[0]= -0.86113631154;
    mGzi4[1]= -0.3399810435;
    mGzi4[2]=  0.3399810435;
    mGzi4[3]=  0.86113631154;
    
    mGzi5[0]= -0.90617984595;
    mGzi5[1]= -0.5384693101;
    mGzi5[2]=  0.0000000000;
    mGzi5[3]=  0.5384693101;
    mGzi5[4]=  0.90617984595;
    
    //積分点での重み 1次元
    mW1[0]= 2.0;

    mW2[0]= 1.0;
    mW2[1]= 1.0;
    
    mW3[0]= 0.5555555555;
    mW3[1]= 0.8888888888;
    mW3[2]= 0.5555555555;
    
    mW4[0]= 0.3478548451;
    mW4[1]= 0.6521451548;
    mW4[2]= 0.6521451548;
    mW4[3]= 0.3478548451;
    
    mW5[0]= 0.2369268851;
    mW5[1]= 0.4786286705;
    mW5[2]= 0.5688888889;
    mW5[3]= 0.4786286705;
    mW5[4]= 0.2369268851;
    
    //積分点数一覧 積分点の個数
    mvIntegNum.resize(5);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 8;
    mvIntegNum[2]= 27;
    mvIntegNum[3]= 64;
    mvIntegNum[4]= 125;
    //積分点数一覧 1D :1次元方向の積分点の個数
    mvIntegNum1D.resize(5);
    mvIntegNum1D[0]= 1;
    mvIntegNum1D[1]= 2;
    mvIntegNum1D[2]= 3;
    mvIntegNum1D[3]= 4;
    mvIntegNum1D[4]= 5;
    
    //--
    //積分点での重み 3次元: w[ir]*w[is]*w[it]
    //--
    mW3d1[0]= mW1[0] * mW1[0] * mW1[0];// 3D:1点 積分の重み => 1点
    setupWeight3d(mW2, 2, mW3d2);// 3D:2点 積分の重み => 8点
    setupWeight3d(mW3, 3, mW3d3);// 3D:3点 積分の重み => 27点
    setupWeight3d(mW4, 4, mW3d4);// 3D:4点 積分の重み => 64点
    setupWeight3d(mW5, 5, mW3d5);// 3D:5点 積分の重み => 125点

    // --
    // 形状関数:
    //  N[積分点][i]
    // --
    // 8節点
    setupShapeFunction(mvN81, 1, 8, mGzi1);//積分点1点(1)
    setupShapeFunction(mvN82, 2, 8, mGzi2);//積分点2点(8)
    // 20節点
    setupShapeFunction(mvN201, 1, 20, mGzi1);// 20節点要素 積分点1
    setupShapeFunction(mvN202, 2, 20, mGzi2);// 20節点要素 積分点2(8)
    setupShapeFunction(mvN203, 3, 20, mGzi3);// 20節点要素 積分点3(27)
//    setupShapeFunction(mvN204, 4, 20, mGzi4);// 20節点要素 積分点4(64)
//    setupShapeFunction(mvN205, 5, 20, mGzi5);// 20節点要素 積分点5(125)
    
    //--
    //自然座標の形状導関数 : dNdr[積分点][i][rst]
    //--
    // 8節点
    setupShapeDeriv(mvdNdr81, 1, 8, mGzi1);// 積分点1
    setupShapeDeriv(mvdNdr82, 2, 8, mGzi2);// 積分点2(8)
    // 20節点
    setupShapeDeriv(mvdNdr201, 1, 20, mGzi1);// 20節点要素 積分点1
    setupShapeDeriv(mvdNdr202, 2, 20, mGzi2);// 20節点要素 積分点2(8)
    setupShapeDeriv(mvdNdr203, 3, 20, mGzi3);// 20節点要素 積分点3(27)
//    setupShapeDeriv(mvdNdr204, 4, 20, mGzi4);// 20節点要素 積分点4(64)
//    setupShapeDeriv(mvdNdr205, 5, 20, mGzi5);// 20節点要素 積分点5(125)


    //--
    //形状関数積分値のセットアップ
    //--
    mvIntegValue20.resize(20);
    mvIntegValue8.resize(8);
    
    setupIntegValue20();
    setupIntegValue8();
}

CShapeHexa::~CShapeHexa()
{
    ;
}


// 1D Weightの3D化
//
// 積分点での重み 3次元: w[ir]*w[is]*w[it]
//
void CShapeHexa::setupWeight3d(const double w[], const uint& numOfIntg, double w3d[])
{
    uint ir,is,it;//積分点座標,各辺(r,s,t)での位置
    uint igauss;  //積分点の連番Index
    //  numOfIntg == 1次元での積分点数

    igauss= 0;
    for(ir=0; ir < numOfIntg; ir++){
    for(is=0; is < numOfIntg; is++){
    for(it=0; it < numOfIntg; it++){
        w3d[igauss]= w[ir]*w[is]*w[it];
        igauss++;
    };
    };
    };
}



//形状関数のセットアップ
//
void CShapeHexa::setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, double Gzi[])
{
    uint igauss;
    uint ir,is,it;
    double r,s,t;
    double r2,s2,t2;

    switch(numOfShape){
        case(8):
            igauss= 0;
            for(ir=0; ir < numOfIntg; ir++){
            for(is=0; is < numOfIntg; is++){
            for(it=0; it < numOfIntg; it++){
                r= Gzi[ir]; s= Gzi[is]; t= Gzi[it];

                ShapeFunction8(N, igauss, r,s,t);
                igauss++;
            };
            };
            };
            break;
        case(20):
            igauss= 0;
            for(ir=0; ir < numOfIntg; ir++){
            for(is=0; is < numOfIntg; is++){
            for(it=0; it < numOfIntg; it++){
                r= Gzi[ir]; s= Gzi[is]; t= Gzi[it];
                r2= r*r;      s2= s*s;      t2= t*t;

                ShapeFunction20(N, igauss, r,s,t, r2,s2,t2);
                igauss++;
            };
            };
            };
            break;
    }
}

//導関数のセットアップ
//
void CShapeHexa::setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, double Gzi[])
{
    uint igauss;
    uint ir,is,it;
    double r,s,t;
    double r2,s2,t2;
    
    switch(numOfShape){
        case(8):
            igauss= 0;
            for(ir=0; ir < numOfIntg; ir++){
            for(is=0; is < numOfIntg; is++){
            for(it=0; it < numOfIntg; it++){
                r= Gzi[ir]; s= Gzi[is]; t= Gzi[it];

                ShapeDeriv8(dNdr, igauss, r,s,t);
                igauss++;
            };
            };
            };
            break;
        case(20):
            igauss= 0;
            for(ir=0; ir < numOfIntg; ir++){
            for(is=0; is < numOfIntg; is++){
            for(it=0; it < numOfIntg; it++){
                r= Gzi[ir]; s= Gzi[is]; t= Gzi[it];
                r2= r*r;      s2= s*s;      t2= t*t;

                ShapeDeriv20(dNdr, igauss, r,s,t, r2,s2,t2);
                igauss++;
            };
            };
            };
            break;
    }
}

// 8 節点の形状関数
//
void CShapeHexa::ShapeFunction8(vvdouble& N, const uint& igauss, 
        const double& r, const double& s, const double& t)
{
    N[igauss][0]= 0.125 * (1.0 - r) * (1.0 - s) * (1.0 - t);
    N[igauss][1]= 0.125 * (1.0 + r) * (1.0 - s) * (1.0 - t);
    N[igauss][2]= 0.125 * (1.0 + r) * (1.0 + s) * (1.0 - t);
    N[igauss][3]= 0.125 * (1.0 - r) * (1.0 + s) * (1.0 - t);
    N[igauss][4]= 0.125 * (1.0 - r) * (1.0 - s) * (1.0 + t);
    N[igauss][5]= 0.125 * (1.0 + r) * (1.0 - s) * (1.0 + t);
    N[igauss][6]= 0.125 * (1.0 + r) * (1.0 + s) * (1.0 + t);
    N[igauss][7]= 0.125 * (1.0 - r) * (1.0 + s) * (1.0 + t);
}
// 20節点の形状関数
//
void CShapeHexa::ShapeFunction20(vvdouble& N,
        const uint& igauss,
        const double& r, const double& s, const double& t,
        const double& r2, const double& s2, const double& t2)
{
    N[igauss][0]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 - t)*(2.0+r+s+t);
    N[igauss][1]= -0.125*(1.0 + r)*(1.0 - s)*(1.0 - t)*(2.0-r+s+t);
    N[igauss][2]= -0.125*(1.0 + r)*(1.0 + s)*(1.0 - t)*(2.0-r-s+t);
    N[igauss][3]= -0.125*(1.0 - r)*(1.0 + s)*(1.0 - t)*(2.0+r-s+t);
    N[igauss][4]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 + t)*(2.0+r+s-t);
    N[igauss][5]= -0.125*(1.0 + r)*(1.0 - s)*(1.0 + t)*(2.0-r+s-t);
    N[igauss][6]= -0.125*(1.0 + r)*(1.0 + s)*(1.0 + t)*(2.0-r-s-t);
    N[igauss][7]= -0.125*(1.0 - r)*(1.0 + s)*(1.0 + t)*(2.0+r-s-t);

    N[igauss][8] = 0.25*(1.0-r2)*(1.0 - s)*(1.0 - t);
    N[igauss][9] = 0.25*(1.0 + r)*(1.0-s2)*(1.0 - t);
    N[igauss][10]= 0.25*(1.0-r2)*(1.0 + s)*(1.0 - t);
    N[igauss][11]= 0.25*(1.0 - r)*(1.0-s2)*(1.0 - t);
    N[igauss][12]= 0.25*(1.0-r2)*(1.0 - s)*(1.0 + t);
    N[igauss][13]= 0.25*(1.0 + r)*(1.0-s2)*(1.0 + t);
    N[igauss][14]= 0.25*(1.0-r2)*(1.0 + s)*(1.0 + t);
    N[igauss][15]= 0.25*(1.0 - r)*(1.0-s2)*(1.0 + t);
    N[igauss][16]= 0.25*(1.0 - r)*(1.0 - s)*(1.0-t2);
    N[igauss][17]= 0.25*(1.0 + r)*(1.0 - s)*(1.0-t2);
    N[igauss][18]= 0.25*(1.0 + r)*(1.0 + s)*(1.0-t2);
    N[igauss][19]= 0.25*(1.0 - r)*(1.0 + s)*(1.0-t2);
}


// 8 節点の形状導関数
void CShapeHexa::ShapeDeriv8(vvvdouble& dNdr, const uint& igauss, 
        const double& r, const double& s, const double& t)
{
    // dN/dr , dN/ds , dN/dt
    // ---------------------
    dNdr[igauss][0][0]= -0.125*(1.0 - s)*(1.0 - t); dNdr[igauss][0][1]= -0.125*(1.0 - r)*(1.0 - t); dNdr[igauss][0][2]= -0.125*(1.0 - r)*(1.0 - s);
    dNdr[igauss][1][0]=  0.125*(1.0 - s)*(1.0 - t); dNdr[igauss][1][1]= -0.125*(1.0 + r)*(1.0 - t); dNdr[igauss][1][2]= -0.125*(1.0 + r)*(1.0 - s);
    dNdr[igauss][2][0]=  0.125*(1.0 + s)*(1.0 - t); dNdr[igauss][2][1]=  0.125*(1.0 + r)*(1.0 - t); dNdr[igauss][2][2]= -0.125*(1.0 + r)*(1.0 + s);
    dNdr[igauss][3][0]= -0.125*(1.0 + s)*(1.0 - t); dNdr[igauss][3][1]=  0.125*(1.0 - r)*(1.0 - t); dNdr[igauss][3][2]= -0.125*(1.0 - r)*(1.0 + s);
    dNdr[igauss][4][0]= -0.125*(1.0 - s)*(1.0 + t); dNdr[igauss][4][1]= -0.125*(1.0 - r)*(1.0 + t); dNdr[igauss][4][2]=  0.125*(1.0 - r)*(1.0 - s);
    dNdr[igauss][5][0]=  0.125*(1.0 - s)*(1.0 + t); dNdr[igauss][5][1]= -0.125*(1.0 + r)*(1.0 + t); dNdr[igauss][5][2]=  0.125*(1.0 + r)*(1.0 - s);
    dNdr[igauss][6][0]=  0.125*(1.0 + s)*(1.0 + t); dNdr[igauss][6][1]=  0.125*(1.0 + r)*(1.0 + t); dNdr[igauss][6][2]=  0.125*(1.0 + r)*(1.0 + s);
    dNdr[igauss][7][0]= -0.125*(1.0 + s)*(1.0 + t); dNdr[igauss][7][1]=  0.125*(1.0 - r)*(1.0 + t); dNdr[igauss][7][2]=  0.125*(1.0 - r)*(1.0 + s);
}

// 20節点の形状導関数
//
void CShapeHexa::ShapeDeriv20(vvvdouble& dNdr,
        const uint& igauss,
        const double& r, const double& s, const double& t,
        const double& r2, const double& s2, const double& t2)
{
    //  dN/dr
    dNdr[igauss][0][0]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 - t)+0.125*(1.0 - s)*(1.0 - t)*(2.0+r+s+t);
    dNdr[igauss][1][0]= +0.125*(1.0 + r)*(1.0 - s)*(1.0 - t)-0.125*(1.0 - s)*(1.0 - t)*(2.0-r+s+t);
    dNdr[igauss][2][0]= +0.125*(1.0 + r)*(1.0 + s)*(1.0 - t)-0.125*(1.0 + s)*(1.0 - t)*(2.0-r-s+t);
    dNdr[igauss][3][0]= -0.125*(1.0 - r)*(1.0 + s)*(1.0 - t)+0.125*(1.0 + s)*(1.0 - t)*(2.0+r-s+t);
    dNdr[igauss][4][0]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 + t)+0.125*(1.0 - s)*(1.0 + t)*(2.0+r+s-t);
    dNdr[igauss][5][0]= +0.125*(1.0 + r)*(1.0 - s)*(1.0 + t)-0.125*(1.0 - s)*(1.0 + t)*(2.0-r+s-t);
    dNdr[igauss][6][0]= +0.125*(1.0 + r)*(1.0 + s)*(1.0 + t)-0.125*(1.0 + s)*(1.0 + t)*(2.0-r-s-t);
    dNdr[igauss][7][0]= -0.125*(1.0 - r)*(1.0 + s)*(1.0 + t)+0.125*(1.0 + s)*(1.0 + t)*(2.0+r-s-t);
    dNdr[igauss][8][0]= -0.50*r*(1.0 - s)*(1.0 - t);
    dNdr[igauss][9][0]= +0.25*(1.0 - s2)*(1.0 - t);
    dNdr[igauss][10][0]= -0.50*r*(1.0 + s)*(1.0 - t);
    dNdr[igauss][11][0]= -0.25*(1.0 - s2)*(1.0 - t);
    dNdr[igauss][12][0]= -0.50*r*(1.0 - s)*(1.0 + t);
    dNdr[igauss][13][0]= +0.25*(1.0 - s2)*(1.0 + t);
    dNdr[igauss][14][0]= -0.50*r*(1.0 + s)*(1.0 + t);
    dNdr[igauss][15][0]= -0.25*(1.0 - s2)*(1.0 + t);
    dNdr[igauss][16][0]= -0.25*(1.0 - s)*(1.0 - t2);
    dNdr[igauss][17][0]= +0.25*(1.0 - s)*(1.0 - t2);
    dNdr[igauss][18][0]= +0.25*(1.0 + s)*(1.0 - t2);
    dNdr[igauss][19][0]= -0.25*(1.0 + s)*(1.0 - t2);

    //  dN/ds
    dNdr[igauss][0][1]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 - t)+0.125*(1.0 - r)*(1.0 - t)*(2.0+r+s+t);
    dNdr[igauss][1][1]= -0.125*(1.0 + r)*(1.0 - s)*(1.0 - t)+0.125*(1.0 + r)*(1.0 - t)*(2.0-r+s+t);
    dNdr[igauss][2][1]= +0.125*(1.0 + r)*(1.0 + s)*(1.0 - t)-0.125*(1.0 + r)*(1.0 - t)*(2.0-r-s+t);
    dNdr[igauss][3][1]= +0.125*(1.0 - r)*(1.0 + s)*(1.0 - t)-0.125*(1.0 - r)*(1.0 - t)*(2.0+r-s+t);
    dNdr[igauss][4][1]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 + t)+0.125*(1.0 - r)*(1.0 + t)*(2.0+r+s-t);
    dNdr[igauss][5][1]= -0.125*(1.0 + r)*(1.0 - s)*(1.0 + t)+0.125*(1.0 + r)*(1.0 + t)*(2.0-r+s-t);
    dNdr[igauss][6][1]= +0.125*(1.0 + r)*(1.0 + s)*(1.0 + t)-0.125*(1.0 + r)*(1.0 + t)*(2.0-r-s-t);
    dNdr[igauss][7][1]= +0.125*(1.0 - r)*(1.0 + s)*(1.0 + t)-0.125*(1.0 - r)*(1.0 + t)*(2.0+r-s-t);
    dNdr[igauss][8][1]= -0.25*(1.0 - r2)*(1.0 - t);
    dNdr[igauss][9][1]= -0.50*(1.0 + r)*s*(1.0 - t);
    dNdr[igauss][10][1]= +0.25*(1.0 - r2)*(1.0 - t);
    dNdr[igauss][11][1]= -0.50*(1.0 - r)*s*(1.0 - t);
    dNdr[igauss][12][1]= -0.25*(1.0 - r2)*(1.0 + t);
    dNdr[igauss][13][1]= -0.50*(1.0 + r)*s*(1.0 + t);
    dNdr[igauss][14][1]= +0.25*(1.0 - r2)*(1.0 + t);
    dNdr[igauss][15][1]= -0.50*(1.0 - r)*s*(1.0 + t);
    dNdr[igauss][16][1]= -0.25*(1.0 - r)*(1.0 - t2);
    dNdr[igauss][17][1]= -0.25*(1.0 + r)*(1.0 - t2);
    dNdr[igauss][18][1]= +0.25*(1.0 + r)*(1.0 - t2);
    dNdr[igauss][19][1]= +0.25*(1.0 - r)*(1.0 - t2);

    //  dN/dt
    dNdr[igauss][0][2]= -0.125*(1.0 - r)*(1.0 - s)*(1.0 - t)+0.125*(1.0 - r)*(1.0 - s)*(2.0+r+s+t);
    dNdr[igauss][1][2]= -0.125*(1.0 + r)*(1.0 - s)*(1.0 - t)+0.125*(1.0 + r)*(1.0 - s)*(2.0-r+s+t);
    dNdr[igauss][2][2]= -0.125*(1.0 + r)*(1.0 + s)*(1.0 - t)+0.125*(1.0 + r)*(1.0 + s)*(2.0-r-s+t);
    dNdr[igauss][3][2]= -0.125*(1.0 - r)*(1.0 + s)*(1.0 - t)+0.125*(1.0 - r)*(1.0 + s)*(2.0+r-s+t);
    dNdr[igauss][4][2]= +0.125*(1.0 - r)*(1.0 - s)*(1.0 + t)-0.125*(1.0 - r)*(1.0 - s)*(2.0+r+s-t);
    dNdr[igauss][5][2]= +0.125*(1.0 + r)*(1.0 - s)*(1.0 + t)-0.125*(1.0 + r)*(1.0 - s)*(2.0-r+s-t);
    dNdr[igauss][6][2]= +0.125*(1.0 + r)*(1.0 + s)*(1.0 + t)-0.125*(1.0 + r)*(1.0 + s)*(2.0-r-s-t);
    dNdr[igauss][7][2]= +0.125*(1.0 - r)*(1.0 + s)*(1.0 + t)-0.125*(1.0 - r)*(1.0 + s)*(2.0+r-s-t);
    dNdr[igauss][8][2]= -0.25*(1.0 - r2)*(1.0 - s);
    dNdr[igauss][9][2]= -0.25*(1.0 + r)*(1.0 - s2);
    dNdr[igauss][10][2]= -0.25*(1.0 - r2)*(1.0 + s);
    dNdr[igauss][11][2]= -0.25*(1.0 - r)*(1.0 - s2);
    dNdr[igauss][12][2]= +0.25*(1.0 - r2)*(1.0 - s);
    dNdr[igauss][13][2]= +0.25*(1.0 + r)*(1.0 - s2);
    dNdr[igauss][14][2]= +0.25*(1.0 - r2)*(1.0 + s);
    dNdr[igauss][15][2]= +0.25*(1.0 - r)*(1.0 - s2);
    dNdr[igauss][16][2]= -0.5*(1.0 - r)*(1.0 - s)*t;
    dNdr[igauss][17][2]= -0.5*(1.0 + r)*(1.0 - s)*t;
    dNdr[igauss][18][2]= -0.5*(1.0 + r)*(1.0 + s)*t;
    dNdr[igauss][19][2]= -0.5*(1.0 - r)*(1.0 + s)*t;
}

//--
//形状関数   引数:積分点位置3次元上のIndex, 形状関数番号(節点)
//--
double& CShapeHexa::N81(const uint& igauss, const uint& ishape){  return mvN81[igauss][ishape];}
double& CShapeHexa::N82(const uint& igauss, const uint& ishape){  return mvN82[igauss][ishape];}
double& CShapeHexa::N201(const uint& igauss, const uint& ishape){ return mvN201[igauss][ishape];}
double& CShapeHexa::N202(const uint& igauss, const uint& ishape){ return mvN202[igauss][ishape];}
double& CShapeHexa::N203(const uint& igauss, const uint& ishape){ return mvN203[igauss][ishape];}
//double& CShapeHexa::N204(const uint& igauss, const uint& ishape){ return mvN204[igauss][ishape];}
//double& CShapeHexa::N205(const uint& igauss, const uint& ishape){ return mvN205[igauss][ishape];}
vdouble& CShapeHexa::N81(const uint& igauss){  return mvN81[igauss];}
vdouble& CShapeHexa::N82(const uint& igauss){  return mvN82[igauss];}
vdouble& CShapeHexa::N201(const uint& igauss){ return mvN201[igauss];}
vdouble& CShapeHexa::N202(const uint& igauss){ return mvN202[igauss];}
vdouble& CShapeHexa::N203(const uint& igauss){ return mvN203[igauss];}

//--
//導関数   引数:積分点位置3次元上のIndex, 形状関数番号(節点),微分方向(axis:0,1,2)
//--
double& CShapeHexa::dNdr81(const uint& igauss, const uint& ishape, const uint& axis){  return mvdNdr81[igauss][ishape][axis];}
double& CShapeHexa::dNdr82(const uint& igauss, const uint& ishape, const uint& axis){  return mvdNdr82[igauss][ishape][axis];}
double& CShapeHexa::dNdr201(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr201[igauss][ishape][axis];}
double& CShapeHexa::dNdr202(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr202[igauss][ishape][axis];}
double& CShapeHexa::dNdr203(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr203[igauss][ishape][axis];}
//double& CShapeHexa::dNdr204(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr204[igauss][ishape][axis];}
//double& CShapeHexa::dNdr205(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdr205[igauss][ishape][axis];}
vvdouble& CShapeHexa::dNdr81(const uint& igauss){ return mvdNdr81[igauss];}
vvdouble& CShapeHexa::dNdr82(const uint& igauss){ return mvdNdr82[igauss];}
vvdouble& CShapeHexa::dNdr201(const uint& igauss){ return mvdNdr201[igauss];}
vvdouble& CShapeHexa::dNdr202(const uint& igauss){ return mvdNdr202[igauss];}
vvdouble& CShapeHexa::dNdr203(const uint& igauss){ return mvdNdr203[igauss];}

// Weight 3D
// --
// 重み配列を返す, 引数:一次元での積分点数
// --
double* CShapeHexa::Weight3d(const uint& integNum1D)
{
    switch(integNum1D){
        case(1):
            return mW3d1;
            break;
        case(2):
            return mW3d2;
            break;
        case(3):
            return mW3d3;
            break;
        case(4):
            return mW3d4;
            break;
        case(5):
            return mW3d5;
            break;
        default:
            break;
    }
}
// Weight 3D
// --
// 重みを直接返す, 引数:配列Index
// --
double& CShapeHexa::Weight3dpt1(){ return mW3d1[0];}
double& CShapeHexa::Weight3dpt2(const uint& igauss){ return mW3d2[igauss];}
double& CShapeHexa::Weight3dpt3(const uint& igauss){ return mW3d3[igauss];}
double& CShapeHexa::Weight3dpt4(const uint& igauss){ return mW3d4[igauss];}
double& CShapeHexa::Weight3dpt5(const uint& igauss){ return mW3d5[igauss];}

// CShapeFunctionBase に移動
//
////  dNdxの計算 3次元(3*3)
////
//void CShapeHexa::ShapeDerivXYZ(const uint& numOfInteg, const uint& numOfShape, const uint& dof,
//                    const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx)
//{
//
//    uint igauss,ishape,row,col;
//    double val;//行列の加算結果 作業変数
//
//    pJacobi->Calculate_J_invJ(dNdr, vNode);// ヤコビアンJ,invJ の計算
//
//    for(igauss=0; igauss< numOfInteg; igauss++){
//    for(ishape=0; ishape< numOfShape; ishape++){
//        // dNdx,dNdy,dNdzの計算
//        for(row=0; row< 3; row++){
//        val=0.0;
//            for(col=0; col< 3; col++){
//                val += dNdr[igauss][ishape][col] * pJacobi->inverse_J(igauss,ishape,row,col);
//            };
//        dNdx[igauss][ishape][row]= val;
//        };
//    };
//    };
//}

// 8節点:Hexa dNdxのセットアップ
// --
void CShapeHexa::Calc_dNdx8(const uint& numOfInteg, CElement *pElement)
//void CShapeHexa::Calc_dNdx8(const uint& numOfInteg, vector<CNode*>& vNode)
{
    vector<CNode*> vNode= pElement->getNode();
    uint numOfShape(8);
    uint dof(3);//3次元: J(3*3)行列

    switch(numOfInteg){
        case(1):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr81, mpJacobi81, vNode, mvdNdx81, mv_detJ81);
            break;
        case(8):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr82, mpJacobi82, vNode, mvdNdx82, mv_detJ82);
            break;
    }
}
// 20節点:Hexa2  dNdxのセットアップ
// --
void CShapeHexa::Calc_dNdx20(const uint& numOfInteg, CElement *pElement)
//void CShapeHexa::Calc_dNdx20(const uint& numOfInteg, vector<CNode*>& vVertNode, vector<CNode*>& vEdgeNode)
{
    vector<CNode*> vNode= pElement->getNode();
    uint numOfShape(20);


//    CNode *pNode;
//    uint ishape,iedge;
//    // 要素の辺ノードをFrontISTR順に並び替えて,vNodeに追加
//    //
//    for(ishape=8; ishape< numOfShape; ishape++){
//        iedge= mpISTR2Edge->HexaEdgeNum(ishape);
//        pNode= pElement->getEdgeInterNode(iedge);
//        vNode.push_back(pNode);
//    };


    uint dof(3);//3次元
    switch(numOfInteg){
        case(1):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr201, mpJacobi201, vNode, mvdNdx201, mv_detJ201);
            break;
        case(8):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr202, mpJacobi202, vNode, mvdNdx202, mv_detJ202);
            break;
        case(27):
            ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr203, mpJacobi203, vNode, mvdNdx203, mv_detJ203);
            break;
    }
}

// detJ の提供
//
double& CShapeHexa::detJ(const uint& elemType, const uint& numOfInteg, const uint& igauss, const uint& ishape)
{
    switch(elemType){
        case(ElementType::Hexa):
            if(numOfInteg==1)
                return mv_detJ81[igauss];
            if(numOfInteg==8)
                return mv_detJ82[igauss];
            break;
        case(ElementType::Hexa2):
            if(numOfInteg==1)
                return mv_detJ201[igauss];
            if(numOfInteg==8)
                return mv_detJ202[igauss];
            if(numOfInteg==27)
                return mv_detJ203[igauss];
            break;
    }
}



// Equivalent Node Force(等価節点力)のための形状関数積分:コンストラクターからコール
// --
// 20節点は,27点積分を使用
// --
void CShapeHexa::setupIntegValue20()
{
    uint ishape;
    for(ishape=0; ishape < 20; ishape++) mvIntegValue20[ishape]= 0.0;

    uint igauss;
    for(igauss=0; igauss < 27; igauss++){
        for(ishape=0; ishape < 20; ishape++){
            mvIntegValue20[ishape] += mvN203[igauss][ishape];
        };
    };

    // Normalize
    for(ishape=0; ishape < 20; ishape++) mvIntegValue20[ishape] /= 27.0;
}
// --
// 8節点
// --
void CShapeHexa::setupIntegValue8()
{
    uint ishape;
    for(ishape=0; ishape < 8; ishape++){
        mvIntegValue8[ishape]= 0.125;
    };
}



















