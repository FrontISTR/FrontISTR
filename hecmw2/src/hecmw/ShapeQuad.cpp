/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeQuad.cpp
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
#include "ShapeQuad.h"
using namespace pmw;
CShapeQuad::CShapeQuad()
{
    mGzi1[0][0]= 0.0;
    mGzi1[0][1]= 0.0;
    mGzi4[0][0]= -0.5773502691896260;
    mGzi4[0][1]= -0.5773502691896260;
    mGzi4[1][0]=  0.5773502691896260;
    mGzi4[1][1]= -0.5773502691896260;
    mGzi4[2][0]= -0.5773502691896260;
    mGzi4[2][1]=  0.5773502691896260;
    mGzi4[3][0]=  0.5773502691896260;
    mGzi4[3][1]=  0.5773502691896260;
    mGzi9[0][0]= -0.7745966692414830;
    mGzi9[0][1]= -0.7745966692414830;
    mGzi9[1][0]= -0.0;
    mGzi9[1][1]= -0.7745966692414830;
    mGzi9[2][0]= 0.7745966692414830;
    mGzi9[2][1]= -0.7745966692414830;
    mGzi9[3][0]= -0.7745966692414830;
    mGzi9[3][1]=  0.0;
    mGzi9[4][0]= -0.0;
    mGzi9[4][1]=  0.0;
    mGzi9[5][0]= 0.7745966692414830;
    mGzi9[5][1]=  0.0;
    mGzi9[6][0]= -0.7745966692414830;
    mGzi9[6][1]=  0.7745966692414830;
    mGzi9[7][0]= -0.0;
    mGzi9[7][1]=  0.7745966692414830;
    mGzi9[8][0]= 0.7745966692414830;
    mGzi9[8][1]=  0.7745966692414830;
    mW1[0]= 4.0;
    mW4[0]= 1.0;
    mW4[1]= 1.0;
    mW4[2]= 1.0;
    mW4[3]= 1.0;
    mW9[0]= 0.3086419753086420;
    mW9[1]= 0.4938271604938270;
    mW9[2]= 0.3086419753086420;
    mW9[3]= 0.4938271604938270;
    mW9[4]= 0.7901234567901230;
    mW9[5]= 0.4938271604938270;
    mW9[6]= 0.3086419753086420;
    mW9[7]= 0.4938271604938270;
    mW9[8]= 0.3086419753086420;
    mvIntegNum.resize(3);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 4;
    mvIntegNum[2]= 9;
    mvIntegNum1D.resize(3);
    mvIntegNum1D[0]= 1;
    mvIntegNum1D[1]= 2;
    mvIntegNum1D[2]= 3;
    uiint numOfIntg,numOfShape,dof;
    numOfIntg=1;
    numOfShape=4;
    dof=2;
    ResizeShape(mvN41, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr41, numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr41, numOfIntg, numOfShape,dof);
    numOfIntg=4;
    numOfShape=8;
    dof=2;
    ResizeShape(mvN84, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr84, numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr84, numOfIntg, numOfShape, dof);
    numOfIntg=9;
    numOfShape=8;
    dof=2;
    ResizeShape(mvN89, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr89, numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr89, numOfIntg, numOfShape, dof);
    setupShapeFunction(mvN41, 1, 4, mGzi1);
    setupShapeFunction(mvN84, 4, 8, mGzi4);
    setupShapeFunction(mvN89, 9, 8, mGzi9);
    setupShapeDeriv(mvdNdr41, 1, 4, mGzi1);
    setupShapeDeriv(mvdNdr84, 4, 8, mGzi4);
    setupShapeDeriv(mvdNdr89, 9, 8, mGzi9);
    setupShape2ndDeriv(mvd2Ndr41, 1, 4, mGzi1);
    setupShape2ndDeriv(mvd2Ndr84, 4, 8, mGzi4);
    setupShape2ndDeriv(mvd2Ndr89, 9, 8, mGzi9);
    mvIntegValue8.resize(8);
    mvIntegValue4.resize(4);
    setupIntegValue8();
    setupIntegValue4();
}
CShapeQuad::~CShapeQuad()
{
    ;
}
void CShapeQuad::setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2])
{
    uiint igauss;
    double r,s;
    switch(numOfShape) {
    case(4):
        igauss= 0;
        r= Gzi[0][0];
        s= Gzi[0][1];
        ShapeFunction4(N, igauss, r,s);
        break;
    case(8):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            r= Gzi[igauss][0];
            s= Gzi[igauss][1];
            ShapeFunction8(N, igauss, r,s);
        };
        break;
    default:
        break;
    }
}
void CShapeQuad::setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2])
{
    uiint igauss;
    double r,s;
    switch(numOfShape) {
    case(4):
        igauss= 0;
        r= Gzi[0][0];
        s= Gzi[0][1];
        ShapeDeriv4(dNdr, igauss, r,s);
        break;
    case(8):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            r= Gzi[igauss][0];
            s= Gzi[igauss][1];
            ShapeDeriv8(dNdr, igauss, r,s);
        };
        break;
    default:
        break;
    }
}
void CShapeQuad::setupShape2ndDeriv(v4double& d2Ndr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2])
{
    uiint igauss;
    double r,s;
    switch(numOfShape) {
    case(4):
        Shape_2ndDeriv4();
        break;
    case(8):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            r= Gzi[igauss][0];
            s= Gzi[igauss][1];
            Shape_2ndDeriv8(d2Ndr, igauss, r,s);
        };
        break;
    default:
        break;
    }
}
void CShapeQuad::ShapeFunction4(vvdouble& N, const uiint& igauss,
                                const double& r, const double& s)
{
    N[igauss][0]= 0.25*(1.0-r)*(1.0-s);
    N[igauss][1]= 0.25*(1.0+r)*(1.0-s);
    N[igauss][2]= 0.25*(1.0+r)*(1.0+s);
    N[igauss][3]= 0.25*(1.0-r)*(1.0+s);
}
void CShapeQuad::ShapeFunction8(vvdouble& N, const uiint& igauss,
                                const double& r, const double& s)
{
    N[igauss][0]= 0.25*(1.0 - r)*(1.0 - s)*(-1.0-r-s);
    N[igauss][1]= 0.25*(1.0 + r)*(1.0 - s)*(-1.0+r-s);
    N[igauss][2]= 0.25*(1.0 + r)*(1.0 + s)*(-1.0+r+s);
    N[igauss][3]= 0.25*(1.0 - r)*(1.0 + s)*(-1.0-r+s);
    N[igauss][4]= 0.50*(1.0-r*r)*(1.0-s);
    N[igauss][5]= 0.50*(1.0-s*s)*(1.0+r);
    N[igauss][6]= 0.50*(1.0-r*r)*(1.0+s);
    N[igauss][7]= 0.50*(1.0-s*s)*(1.0-r);
}
void CShapeQuad::ShapeDeriv4(vvvdouble& dNdr, const uiint& igauss,
                             const double& r, const double& s)
{
    dNdr[igauss][0][0]=  -0.25*(1.0-s);
    dNdr[igauss][1][0]=   0.25*(1.0-s);
    dNdr[igauss][2][0]=   0.25*(1.0+s);
    dNdr[igauss][3][0]=  -0.25*(1.0+s);
    dNdr[igauss][0][1]=  -0.25*(1.0-r);
    dNdr[igauss][1][1]=  -0.25*(1.0+r);
    dNdr[igauss][2][1]=   0.25*(1.0+r);
    dNdr[igauss][3][1]=   0.25*(1.0-r);
}
void CShapeQuad::ShapeDeriv8(vvvdouble& dNdr, const uiint& igauss,
                             const double& r, const double& s)
{
    dNdr[igauss][0][0]=  0.25*(1.0-s)*(2.0*r+s);
    dNdr[igauss][1][0]=  0.25*(1.0-s)*(2.0*r-s);
    dNdr[igauss][2][0]=  0.25*(1.0+s)*(2.0*r+s);
    dNdr[igauss][3][0]=  0.25*(1.0+s)*(2.0*r-s);
    dNdr[igauss][4][0]= -r*(1.0-s);
    dNdr[igauss][5][0]=  0.50*(1.0-s*s);
    dNdr[igauss][6][0]= -r*(1.0+s);
    dNdr[igauss][7][0]= -0.50*(1.0-s*s);
    dNdr[igauss][0][1]=  0.25*(1.0-r)*(r+2.0*s);
    dNdr[igauss][1][1]= -0.25*(1.0+r)*(r-2.0*s);
    dNdr[igauss][2][1]=  0.25*(1.0+r)*(r+2.0*s);
    dNdr[igauss][3][1]= -0.25*(1.0-r)*(r-2.0*s);
    dNdr[igauss][4][1]= -0.50*(1.0-r*r);
    dNdr[igauss][5][1]= -s*(1.0+r);
    dNdr[igauss][6][1]=  0.50*(1.0-r*r);
    dNdr[igauss][7][1]= -s*(1.0-r);
}
void CShapeQuad::Shape_2ndDeriv4()
{
    uiint igauss(0);
    uiint ishape;
    for(ishape=0; ishape < 4; ishape++) {
        mvd2Ndr41[igauss][ishape][0][0]= 0.0;
        mvd2Ndr41[igauss][ishape][1][1]= 0.0;
    };
    mvd2Ndr41[igauss][0][0][1]=  0.25;
    mvd2Ndr41[igauss][1][0][1]= -0.25;
    mvd2Ndr41[igauss][2][0][1]=  0.25;
    mvd2Ndr41[igauss][3][0][1]= -0.25;
    mvd2Ndr41[igauss][0][1][0]=  0.25;
    mvd2Ndr41[igauss][1][1][0]= -0.25;
    mvd2Ndr41[igauss][2][1][0]=  0.25;
    mvd2Ndr41[igauss][3][1][0]= -0.25;
}
void CShapeQuad::Shape_2ndDeriv8(v4double& d2Ndr, const uiint& igauss,
                                 const double& r, const double& s)
{
    d2Ndr[igauss][0][0][0]=  0.50*(1.0-s);
    d2Ndr[igauss][1][0][0]=  0.50*(1.0-s);
    d2Ndr[igauss][2][0][0]=  0.50*(1.0+s);
    d2Ndr[igauss][3][0][0]=  0.50*(1.0+s);
    d2Ndr[igauss][4][0][0]=  s-1.0;
    d2Ndr[igauss][5][0][0]=  0.0;
    d2Ndr[igauss][6][0][0]= -s-1.0;
    d2Ndr[igauss][7][0][0]=  0.0;
    d2Ndr[igauss][0][0][1]=  0.25-0.50*(r+s);
    d2Ndr[igauss][1][0][1]= -0.25-0.50*(r-s);
    d2Ndr[igauss][2][0][1]=  0.25+0.50*(r+s);
    d2Ndr[igauss][3][0][1]= -0.25+0.50*(r-s);
    d2Ndr[igauss][4][0][1]=  r;
    d2Ndr[igauss][5][0][1]= -s;
    d2Ndr[igauss][6][0][1]= -r;
    d2Ndr[igauss][7][0][1]=  s;
    d2Ndr[igauss][0][1][0]=  0.25-0.50*(r+s);
    d2Ndr[igauss][1][1][0]= -0.25-0.50*(r-s);
    d2Ndr[igauss][2][1][0]=  0.25+0.50*(r+s);
    d2Ndr[igauss][3][1][0]= -0.25+0.50*(r-s);
    d2Ndr[igauss][4][1][0]=  r;
    d2Ndr[igauss][5][1][0]= -s;
    d2Ndr[igauss][6][1][0]= -r;
    d2Ndr[igauss][7][1][0]=  s;
    d2Ndr[igauss][0][1][1]=  0.50*(1.0-r);
    d2Ndr[igauss][1][1][1]=  0.50*(1.0+r);
    d2Ndr[igauss][2][1][1]=  0.50*(1.0+r);
    d2Ndr[igauss][3][1][1]=  0.50*(1.0-r);
    d2Ndr[igauss][4][1][1]=  0.0;
    d2Ndr[igauss][5][1][1]= -r-1.0;
    d2Ndr[igauss][6][1][1]=  0.0;
    d2Ndr[igauss][7][1][1]=  r-1.0;
}
double& CShapeQuad::N41(const uiint& igauss, const uiint& ishape)
{
    return mvN41[igauss][ishape];
}
double& CShapeQuad::N84(const uiint& igauss, const uiint& ishape)
{
    return mvN84[igauss][ishape];
}
double& CShapeQuad::N89(const uiint& igauss, const uiint& ishape)
{
    return mvN89[igauss][ishape];
}
double& CShapeQuad::dNdr41(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr41[igauss][ishape][axis];
}
double& CShapeQuad::dNdr84(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr84[igauss][ishape][axis];
}
double& CShapeQuad::dNdr89(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr89[igauss][ishape][axis];
}
double& CShapeQuad::d2Ndr41(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1)
{
    return mvd2Ndr41[igauss][ishape][deriv_axis0][deriv_axis1];
}
double& CShapeQuad::d2Ndr84(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1)
{
    return mvd2Ndr84[igauss][ishape][deriv_axis0][deriv_axis1];
}
double& CShapeQuad::d2Ndr89(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1)
{
    return mvd2Ndr89[igauss][ishape][deriv_axis0][deriv_axis1];
}
double* CShapeQuad::Weight(const uiint& integNum)
{
    switch(integNum) {
    case(1):
        return mW1;
        break;
    case(4):
        return mW4;
        break;
    case(9):
        return mW9;
        break;
    default:
        break;
    }
}
double& CShapeQuad::Weight_pt1()
{
    return mW1[0];
}
double& CShapeQuad::Weight_pt4(const uiint& igauss)
{
    return mW4[igauss];
}
double& CShapeQuad::Weight_pt9(const uiint& igauss)
{
    return mW9[igauss];
}
void CShapeQuad::setupIntegValue8()
{
    uiint ishape;
    for(ishape=0; ishape < 4; ishape++) mvIntegValue8[ishape]= -1.0/12.0;
    for(ishape=4; ishape < 8; ishape++) mvIntegValue8[ishape]=  1.0/3.0;
}
void CShapeQuad::setupIntegValue4()
{
    uiint ishape;
    for(ishape=0; ishape < 4; ishape++) mvIntegValue4[ishape]= 0.25;
}
