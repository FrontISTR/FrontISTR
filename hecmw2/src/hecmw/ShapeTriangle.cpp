//
//  ShapeTriangle.cpp
//
//
//
//
//              2010.02.02
//              k.Takeda
#include "ShapeTriangle.h"
using namespace pmw;

CShapeTriangle::CShapeTriangle()
{
    
    
    //積分座標: gzi,eta
    //----
    //1点 積分
    mGzi1[0][0]= 0.333333333333333;  mGzi1[0][1]= 0.333333333333333;
    //3点 積分
    mGzi3[0][0]= 0.166666666666667;  mGzi3[0][1]= 0.166666666666667;
    mGzi3[1][0]= 0.666666666666667;  mGzi3[1][1]= 0.166666666666667;
    mGzi3[2][0]= 0.166666666666667;  mGzi3[2][1]= 0.666666666666667;
    
    //積分点の重み
    //----
    //1点 積分
    mW1[0]= 0.5;
    //3点 積分
    mW3[0]= 0.166666666666666;
    mW3[1]= 0.166666666666666;
    mW3[2]= 0.166666666666666;

    //積分点数
    mvIntegNum.resize(2);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 3;


    uint numOfIntg,numOfShape,dof;
    //領域確保
    numOfIntg=1; numOfShape=3; dof=2;
    ResizeShape(mvN31, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr31,numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr31, numOfIntg, numOfShape, dof);

    numOfIntg=3; numOfShape=6; dof=2;
    ResizeShape(mvN63, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr63,numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr63, numOfIntg, numOfShape, dof);


    //形状関数 & 導関数のセットアップ
    setupShapeFunction();
    setupShapeDeriv();
    setupShape2ndDeriv();


    //形状関数の積分 セットアップ
    mvIntegValue6.resize(6);
    mvIntegValue3.resize(3);
    setupIntegValue6();
    setupIntegValue3();
}
CShapeTriangle::~CShapeTriangle()
{
    ;
}

//形状関数 & 導関数のセットアップ
void CShapeTriangle::setupShapeFunction()
{
    uint igauss;
    uint numOfIntg;
    double r,s;

    igauss= 0;
    r= mGzi1[igauss][0]; s= mGzi1[igauss][1];
    ShapeFunction3(mvN31, igauss, r,s);

    numOfIntg= 3;
    for(igauss=0; igauss < numOfIntg; igauss++){
        r= mGzi3[igauss][0];  s= mGzi3[igauss][1];
        
        ShapeFunction6(mvN63, igauss, r,s);
    };
}
void CShapeTriangle::setupShapeDeriv()
{
    uint igauss;
    uint numOfIntg;
    double r,s;

    igauss= 0;
    ShapeDeriv3(mvdNdr31, igauss);

    numOfIntg= 3;
    for(igauss=0; igauss < numOfIntg; igauss++){
        r= mGzi3[igauss][0];  s= mGzi3[igauss][1];

        ShapeDeriv6(mvdNdr63, igauss, r,s);
    };
}
void CShapeTriangle::setupShape2ndDeriv()
{
    uint igauss;
    uint numOfIntg;

    Shape_2ndDeriv3();

    numOfIntg= 3;
    for(igauss=0; igauss < numOfIntg; igauss++){
        Shape_2ndDeriv6(mvd2Ndr63, igauss);
    };
}


//形状関数
//void CShapeTriangle::ShapeFunction3(double N[][3], const uint& igauss, const double& r, const double& s)
void CShapeTriangle::ShapeFunction3(vvdouble& N, const uint& igauss, const double& r, const double& s)
{
    N[igauss][0]= r;  N[igauss][1]= s;  N[igauss][2]= 1.0 - r - s;
}
//void CShapeTriangle::ShapeFunction6(double N[][6], const uint& igauss, const double& r, const double& s)
void CShapeTriangle::ShapeFunction6(vvdouble& N, const uint& igauss, const double& r, const double& s)
{
    double t= 1.0-r-s;

    N[igauss][0]= t*(2.0*t-1.0);
    N[igauss][1]= r*(2.0*r-1.0);
    N[igauss][2]= s*(2.0*s-1.0);
    N[igauss][3]= 4.0*r*t;
    N[igauss][4]= 4.0*r*s;
    N[igauss][5]= 4.0*s*t;
}

//導関数
//void CShapeTriangle::ShapeDeriv3(double dNdr[][3][2], const uint& igauss)
void CShapeTriangle::ShapeDeriv3(vvvdouble& dNdr, const uint& igauss)
{
    dNdr[igauss][0][0]=  1.0;
    dNdr[igauss][1][0]=  0.0;
    dNdr[igauss][2][0]= -1.0;

    dNdr[igauss][0][1]=  0.0;
    dNdr[igauss][1][1]=  1.0;
    dNdr[igauss][2][1]= -1.0;
}
//void CShapeTriangle::ShapeDeriv6(double dNdr[][6][2], const uint& igauss, const double& r, const double& s)
void CShapeTriangle::ShapeDeriv6(vvvdouble& dNdr, const uint& igauss, const double& r, const double& s)
{
    double t= 1.0-r-s;

    dNdr[igauss][0][0]=  1.0-4.0*t;
    dNdr[igauss][1][0]=  4.0*r-1.0;
    dNdr[igauss][2][0]=  0.0;
    dNdr[igauss][3][0]=  4.0*(t-r);
    dNdr[igauss][4][0]=  4.0*s;
    dNdr[igauss][5][0]= -4.0*s;

    dNdr[igauss][0][1]=  1.0-4.0*t;
    dNdr[igauss][1][1]=  0.0;
    dNdr[igauss][2][1]=  4.0*s-1.0;
    dNdr[igauss][3][1]= -4.0*r;
    dNdr[igauss][4][1]=  4.0*r;
    dNdr[igauss][5][1]=  4.0*(t-s);
}

//2次導関数
void CShapeTriangle::Shape_2ndDeriv3()
{
    uint ishape,iaxis,jaxis;

    for(ishape=0; ishape < 3; ishape++){
    for(iaxis=0;  iaxis  < 2; iaxis++) {
    for(jaxis=0;  jaxis  < 2; jaxis++) {
        mvd2Ndr31[0][ishape][iaxis][jaxis]= 0.0;
    };
    };
    };
}
//void CShapeTriangle::Shape_2ndDeriv6(double d2Ndr[][6][2][2], const uint& igauss)
void CShapeTriangle::Shape_2ndDeriv6(v4double& d2Ndr, const uint& igauss)
{
    d2Ndr[igauss][0][0][0]=  4.0;  d2Ndr[igauss][0][0][1]=  4.0;
    d2Ndr[igauss][1][0][0]=  4.0;  d2Ndr[igauss][1][0][1]=  0.0;
    d2Ndr[igauss][2][0][0]=  0.0;  d2Ndr[igauss][2][0][1]=  0.0;
    d2Ndr[igauss][3][0][0]= -8.0;  d2Ndr[igauss][3][0][1]= -4.0;
    d2Ndr[igauss][4][0][0]=  0.0;  d2Ndr[igauss][4][0][1]=  4.0;
    d2Ndr[igauss][5][0][0]=  0.0;  d2Ndr[igauss][5][0][1]= -4.0;

    d2Ndr[igauss][0][1][0]=  4.0;  d2Ndr[igauss][0][1][1]=  4.0;
    d2Ndr[igauss][1][1][0]=  0.0;  d2Ndr[igauss][1][1][1]=  0.0;
    d2Ndr[igauss][2][1][0]=  0.0;  d2Ndr[igauss][2][1][1]=  4.0;
    d2Ndr[igauss][3][1][0]= -4.0;  d2Ndr[igauss][3][1][1]=  0.0;
    d2Ndr[igauss][4][1][0]=  4.0;  d2Ndr[igauss][4][1][1]=  0.0;
    d2Ndr[igauss][5][1][0]= -4.0;  d2Ndr[igauss][5][1][1]= -8.0;
}


// 形状関数   引数:積分点, 形状関数
// ----
double& CShapeTriangle::N31(const uint& igauss, const uint& ishape){ return mvN31[igauss][ishape];}
double& CShapeTriangle::N63(const uint& igauss, const uint& ishape){ return mvN63[igauss][ishape];}

// 導関数   引数:積分点, 形状関数,座標方向
// ----
double& CShapeTriangle::dNdr31(const uint& igauss, const uint& ishape, const uint& deriv_axis){
    return mvdNdr31[igauss][ishape][deriv_axis];
}
double& CShapeTriangle::dNdr63(const uint& igauss, const uint& ishape, const uint& deriv_axis){
    return mvdNdr63[igauss][ishape][deriv_axis];
}

// 2次導関数 引数:積分点, 形状関数, 座標方向, 座標方向
// ----
double& CShapeTriangle::d2Ndr31(const uint& igauss, const uint& ishape, const uint& iaxis, const uint& jaxis){
    return mvd2Ndr31[igauss][ishape][iaxis][jaxis];
}
double& CShapeTriangle::d2Ndr63(const uint& igauss, const uint& ishape, const uint& iaxis, const uint& jaxis){
    return mvd2Ndr63[igauss][ishape][iaxis][jaxis];
}



// 返り値:重みの配列, 引数:積分点数
double* CShapeTriangle::Weight(const uint& integNum)
{
    switch(integNum){
        case(1):
            return mW1;
            break;
        case(3):
            return mW3;
            break;
        default:
            break;
    }

}
// 返り値:重み, 引数:重みの配列Index
double& CShapeTriangle::Weight_pt1(){ return mW1[0];}
double& CShapeTriangle::Weight_pt3(const uint& igauss){ return mW3[igauss];}



// 形状関数の積分値セットアップ：コンストラクタからコール
// --
// 6節点
// --
void CShapeTriangle::setupIntegValue6()
{
    uint ishape;
    // 頂点
    for(ishape=0; ishape < 3; ishape++) mvIntegValue6[ishape]= 0.0;
    // 辺
    for(ishape=3; ishape < 6; ishape++) mvIntegValue6[ishape]= 1.0/3.0;
}
// --
// 3節点
// --
void CShapeTriangle::setupIntegValue3()
{
    uint ishape;
    for(ishape=0; ishape < 3; ishape++) mvIntegValue3[ishape]= 1.0/3.0;
}







