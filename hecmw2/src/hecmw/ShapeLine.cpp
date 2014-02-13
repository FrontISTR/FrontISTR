/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeLine.cpp
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
#include "ShapeLine.h"
using namespace pmw;
CShapeLine::CShapeLine()
{
    mGzi1[0][0]= 0.0;
    mGzi2[0][0]= -0.577350269189626;
    mGzi2[1][0]=  0.577350269189626;
    mW1[0]= 2.0;
    mW2[0]= 1.0;
    mW2[1]= 1.0;
    mvIntegNum.resize(2);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 2;
    uiint numOfIntg,numOfShape,dof;
    numOfIntg=1;
    numOfShape=2;
    dof=1;
    ResizeShape(mvN21, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr21,numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr21, numOfIntg, numOfShape, dof);
    numOfIntg=2;
    numOfShape=3;
    dof=1;
    ResizeShape(mvN32, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr32,numOfIntg, numOfShape, dof);
    Resize2ndDeriv(mvd2Ndr32, numOfIntg, numOfShape, dof);
    setupShapeFunction();
    setupShapeDeriv();
    setupShape2ndDeriv();
    mvIntegValue3.resize(3);
    mvIntegValue2.resize(2);
    setupIntegValue3();
    setupIntegValue2();
}
CShapeLine::~CShapeLine()
{
    ;
}
void CShapeLine::setupShapeFunction()
{
    uiint igauss;
    double r;
    igauss=0;
    r= mGzi1[0][0];
    ShapeFunction2(mvN21, igauss, r);
    uiint numOfIntg(2);
    for(igauss=0; igauss< numOfIntg; igauss++) {
        r= mGzi2[igauss][0];
        ShapeFunction3(mvN32, igauss, r);
    };
}
void CShapeLine::setupShapeDeriv()
{
    ShapeDeriv2();
    uiint igauss;
    double r;
    uiint numOfIntg(2);
    for(igauss=0; igauss< numOfIntg; igauss++) {
        r= mGzi2[igauss][0];
        ShapeDeriv3(mvdNdr32, igauss, r);
    };
}
void CShapeLine::setupShape2ndDeriv()
{
    Shape_2ndDeriv2();
    Shape_2ndDeriv3();
}
void CShapeLine::ShapeFunction2(vvdouble& N, const uiint& igauss, const double& r)
{
    N[igauss][0]= 0.50*(1.0-r);
    N[igauss][1]= 0.50*(1.0+r);
}
void CShapeLine::ShapeFunction3(vvdouble& N, const uiint& igauss, const double& r)
{
    N[igauss][0]= -0.50*(1.0-r)*r;
    N[igauss][1]=  0.50*(1.0+r)*r;
    N[igauss][2]=  1.00-r*r;
}
void CShapeLine::ShapeDeriv2()
{
    uiint igauss(0);
    mvdNdr21[igauss][0][0]= -0.50;
    mvdNdr21[igauss][1][0]=  0.50;
}
void CShapeLine::ShapeDeriv3(vvvdouble& dNdr, const uiint& igauss, const double& r)
{
    dNdr[igauss][0][0]=  r-0.50;
    dNdr[igauss][1][0]=  r+0.50;
    dNdr[igauss][2][0]= -2.0*r;
}
void CShapeLine::Shape_2ndDeriv2()
{
    uiint igauss(0);
    mvd2Ndr21[igauss][0][0][0]= 0.0;
    mvd2Ndr21[igauss][1][0][0]= 0.0;
}
void CShapeLine::Shape_2ndDeriv3()
{
    uiint igauss(0);
    mvd2Ndr32[igauss][0][0][0]=  1.0;
    mvd2Ndr32[igauss][1][0][0]=  1.0;
    mvd2Ndr32[igauss][2][0][0]= -2.0;
}
double& CShapeLine::N21(const uiint& igauss, const uiint& ishape)
{
    return mvN21[igauss][ishape];
}
double& CShapeLine::N32(const uiint& igauss, const uiint& ishape)
{
    return mvN32[igauss][ishape];
}
double& CShapeLine::dNdr21(const uiint& igauss, const uiint& ishape)
{
    return mvdNdr21[igauss][ishape][0];
}
double& CShapeLine::dNdr32(const uiint& igauss, const uiint& ishape)
{
    return mvdNdr32[igauss][ishape][0];
}
double& CShapeLine::d2Ndr21(const uiint& igauss, const uiint& ishape)
{
    return mvd2Ndr21[igauss][ishape][0][0];
}
double& CShapeLine::d2Ndr32(const uiint& igauss, const uiint& ishape)
{
    return mvd2Ndr32[igauss][ishape][0][0];
}
double* CShapeLine::Weight(const uiint& integNum)
{
    switch(integNum) {
    case(1):
        return mW1;
        break;
    case(2):
        return mW2;
        break;
    default:
        break;
    }
}
double& CShapeLine::Weight_pt1()
{
    return mW1[0];
}
double& CShapeLine::Weight_pt2(const uiint& igauss)
{
    return mW2[igauss];
}
void CShapeLine::setupIntegValue3()
{
    uiint ishape;
    for(ishape=0; ishape < 2; ishape++) mvIntegValue3[ishape]= 1.0/6.0;
    mvIntegValue3[2]= 2.0/3.0;
}
void CShapeLine::setupIntegValue2()
{
    uiint ishape;
    for(ishape=0; ishape < 2; ishape++) mvIntegValue2[ishape]= 0.5;
}
