/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeHexa.cpp
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
#include <math.h>
#include <vector>
#include "ShapeHexa.h"
using namespace pmw;
CShapeHexa::CShapeHexa()
{
    uiint numOfIntg, numOfShape, dof;
    numOfIntg= 1;
    numOfShape= 8;
    dof= 3;
    ResizeShape(mvN81, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr81, numOfIntg, numOfShape, dof);
    numOfIntg= 8;
    numOfShape= 8;
    dof= 3;
    ResizeShape(mvN82, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr82, numOfIntg, numOfShape, dof);
    numOfIntg= 1;
    numOfShape= 20;
    dof= 3;
    ResizeShape(mvN201,  numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr201, numOfIntg, numOfShape, dof);
    numOfIntg= 8;
    numOfShape= 20;
    dof= 3;
    ResizeShape(mvN202,  numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr202, numOfIntg, numOfShape, dof);
    numOfIntg= 27;
    numOfShape= 20;
    dof= 3;
    ResizeShape(mvN203, numOfIntg, numOfShape);
    ResizeDeriv(mvdNdr203, numOfIntg, numOfShape, dof);
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
    numOfIntg= 1;
    numOfShape= 8;
    dof= 3;
    ResizeDeriv(mvdNdx81, numOfIntg, numOfShape, dof);
    numOfIntg= 8;
    numOfShape= 8;
    dof= 3;
    ResizeDeriv(mvdNdx82, numOfIntg, numOfShape, dof);
    numOfIntg= 1;
    numOfShape= 20;
    dof= 3;
    ResizeDeriv(mvdNdx201, numOfIntg, numOfShape, dof);
    numOfIntg= 8;
    numOfShape= 20;
    dof= 3;
    ResizeDeriv(mvdNdx202, numOfIntg, numOfShape, dof);
    numOfIntg= 27;
    numOfShape= 20;
    dof= 3;
    ResizeDeriv(mvdNdx203, numOfIntg, numOfShape, dof);
    numOfIntg= 1;
    numOfShape= 8;
    Resize_detJ(mv_detJ81, numOfIntg);
    numOfIntg= 8;
    numOfShape= 8;
    Resize_detJ(mv_detJ82, numOfIntg);
    numOfIntg= 1;
    numOfShape= 20;
    Resize_detJ(mv_detJ201, numOfIntg);
    numOfIntg= 8;
    numOfShape= 20;
    Resize_detJ(mv_detJ202, numOfIntg);
    numOfIntg= 27;
    numOfShape= 20;
    Resize_detJ(mv_detJ203, numOfIntg);
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
    mvIntegNum.resize(5);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 8;
    mvIntegNum[2]= 27;
    mvIntegNum[3]= 64;
    mvIntegNum[4]= 125;
    mvIntegNum1D.resize(5);
    mvIntegNum1D[0]= 1;
    mvIntegNum1D[1]= 2;
    mvIntegNum1D[2]= 3;
    mvIntegNum1D[3]= 4;
    mvIntegNum1D[4]= 5;
    mW3d1[0]= mW1[0] * mW1[0] * mW1[0];
    setupWeight3d(mW2, 2, mW3d2);
    setupWeight3d(mW3, 3, mW3d3);
    setupWeight3d(mW4, 4, mW3d4);
    setupWeight3d(mW5, 5, mW3d5);
    setupShapeFunction(mvN81, 1, 8, mGzi1);
    setupShapeFunction(mvN82, 2, 8, mGzi2);
    setupShapeFunction(mvN201, 1, 20, mGzi1);
    setupShapeFunction(mvN202, 2, 20, mGzi2);
    setupShapeFunction(mvN203, 3, 20, mGzi3);
    setupShapeDeriv(mvdNdr81, 1, 8, mGzi1);
    setupShapeDeriv(mvdNdr82, 2, 8, mGzi2);
    setupShapeDeriv(mvdNdr201, 1, 20, mGzi1);
    setupShapeDeriv(mvdNdr202, 2, 20, mGzi2);
    setupShapeDeriv(mvdNdr203, 3, 20, mGzi3);
    mvIntegValue20.resize(20);
    mvIntegValue8.resize(8);
    setupIntegValue20();
    setupIntegValue8();
}
CShapeHexa::~CShapeHexa()
{
    ;
}
void CShapeHexa::setupWeight3d(const double w[], const uiint& numOfIntg, double w3d[])
{
    uiint ir,is,it;
    uiint igauss;
    igauss= 0;
    for(ir=0; ir < numOfIntg; ir++) {
        for(is=0; is < numOfIntg; is++) {
            for(it=0; it < numOfIntg; it++) {
                w3d[igauss]= w[ir]*w[is]*w[it];
                igauss++;
            };
        };
    };
}
void CShapeHexa::setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, double Gzi[])
{
    uiint igauss;
    uiint ir,is,it;
    double r,s,t;
    double r2,s2,t2;
    switch(numOfShape) {
    case(8):
        igauss= 0;
        for(ir=0; ir < numOfIntg; ir++) {
            for(is=0; is < numOfIntg; is++) {
                for(it=0; it < numOfIntg; it++) {
                    r= Gzi[ir];
                    s= Gzi[is];
                    t= Gzi[it];
                    ShapeFunction8(N, igauss, r,s,t);
                    igauss++;
                };
            };
        };
        break;
    case(20):
        igauss= 0;
        for(ir=0; ir < numOfIntg; ir++) {
            for(is=0; is < numOfIntg; is++) {
                for(it=0; it < numOfIntg; it++) {
                    r= Gzi[ir];
                    s= Gzi[is];
                    t= Gzi[it];
                    r2= r*r;
                    s2= s*s;
                    t2= t*t;
                    ShapeFunction20(N, igauss, r,s,t, r2,s2,t2);
                    igauss++;
                };
            };
        };
        break;
    }
}
void CShapeHexa::setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, double Gzi[])
{
    uiint igauss;
    uiint ir,is,it;
    double r,s,t;
    double r2,s2,t2;
    switch(numOfShape) {
    case(8):
        igauss= 0;
        for(ir=0; ir < numOfIntg; ir++) {
            for(is=0; is < numOfIntg; is++) {
                for(it=0; it < numOfIntg; it++) {
                    r= Gzi[ir];
                    s= Gzi[is];
                    t= Gzi[it];
                    ShapeDeriv8(dNdr, igauss, r,s,t);
                    igauss++;
                };
            };
        };
        break;
    case(20):
        igauss= 0;
        for(ir=0; ir < numOfIntg; ir++) {
            for(is=0; is < numOfIntg; is++) {
                for(it=0; it < numOfIntg; it++) {
                    r= Gzi[ir];
                    s= Gzi[is];
                    t= Gzi[it];
                    r2= r*r;
                    s2= s*s;
                    t2= t*t;
                    ShapeDeriv20(dNdr, igauss, r,s,t, r2,s2,t2);
                    igauss++;
                };
            };
        };
        break;
    }
}
void CShapeHexa::ShapeFunction8(vvdouble& N, const uiint& igauss,
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
void CShapeHexa::ShapeFunction20(vvdouble& N,
                                 const uiint& igauss,
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
void CShapeHexa::ShapeDeriv8(vvvdouble& dNdr, const uiint& igauss,
                             const double& r, const double& s, const double& t)
{
    dNdr[igauss][0][0]= -0.125*(1.0 - s)*(1.0 - t);
    dNdr[igauss][0][1]= -0.125*(1.0 - r)*(1.0 - t);
    dNdr[igauss][0][2]= -0.125*(1.0 - r)*(1.0 - s);
    dNdr[igauss][1][0]=  0.125*(1.0 - s)*(1.0 - t);
    dNdr[igauss][1][1]= -0.125*(1.0 + r)*(1.0 - t);
    dNdr[igauss][1][2]= -0.125*(1.0 + r)*(1.0 - s);
    dNdr[igauss][2][0]=  0.125*(1.0 + s)*(1.0 - t);
    dNdr[igauss][2][1]=  0.125*(1.0 + r)*(1.0 - t);
    dNdr[igauss][2][2]= -0.125*(1.0 + r)*(1.0 + s);
    dNdr[igauss][3][0]= -0.125*(1.0 + s)*(1.0 - t);
    dNdr[igauss][3][1]=  0.125*(1.0 - r)*(1.0 - t);
    dNdr[igauss][3][2]= -0.125*(1.0 - r)*(1.0 + s);
    dNdr[igauss][4][0]= -0.125*(1.0 - s)*(1.0 + t);
    dNdr[igauss][4][1]= -0.125*(1.0 - r)*(1.0 + t);
    dNdr[igauss][4][2]=  0.125*(1.0 - r)*(1.0 - s);
    dNdr[igauss][5][0]=  0.125*(1.0 - s)*(1.0 + t);
    dNdr[igauss][5][1]= -0.125*(1.0 + r)*(1.0 + t);
    dNdr[igauss][5][2]=  0.125*(1.0 + r)*(1.0 - s);
    dNdr[igauss][6][0]=  0.125*(1.0 + s)*(1.0 + t);
    dNdr[igauss][6][1]=  0.125*(1.0 + r)*(1.0 + t);
    dNdr[igauss][6][2]=  0.125*(1.0 + r)*(1.0 + s);
    dNdr[igauss][7][0]= -0.125*(1.0 + s)*(1.0 + t);
    dNdr[igauss][7][1]=  0.125*(1.0 - r)*(1.0 + t);
    dNdr[igauss][7][2]=  0.125*(1.0 - r)*(1.0 + s);
}
void CShapeHexa::ShapeDeriv20(vvvdouble& dNdr,
                              const uiint& igauss,
                              const double& r, const double& s, const double& t,
                              const double& r2, const double& s2, const double& t2)
{
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
double& CShapeHexa::N81(const uiint& igauss, const uiint& ishape)
{
    return mvN81[igauss][ishape];
}
double& CShapeHexa::N82(const uiint& igauss, const uiint& ishape)
{
    return mvN82[igauss][ishape];
}
double& CShapeHexa::N201(const uiint& igauss, const uiint& ishape)
{
    return mvN201[igauss][ishape];
}
double& CShapeHexa::N202(const uiint& igauss, const uiint& ishape)
{
    return mvN202[igauss][ishape];
}
double& CShapeHexa::N203(const uiint& igauss, const uiint& ishape)
{
    return mvN203[igauss][ishape];
}
vdouble& CShapeHexa::N81(const uiint& igauss)
{
    return mvN81[igauss];
}
vdouble& CShapeHexa::N82(const uiint& igauss)
{
    return mvN82[igauss];
}
vdouble& CShapeHexa::N201(const uiint& igauss)
{
    return mvN201[igauss];
}
vdouble& CShapeHexa::N202(const uiint& igauss)
{
    return mvN202[igauss];
}
vdouble& CShapeHexa::N203(const uiint& igauss)
{
    return mvN203[igauss];
}
double& CShapeHexa::dNdr81(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr81[igauss][ishape][axis];
}
double& CShapeHexa::dNdr82(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr82[igauss][ishape][axis];
}
double& CShapeHexa::dNdr201(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr201[igauss][ishape][axis];
}
double& CShapeHexa::dNdr202(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr202[igauss][ishape][axis];
}
double& CShapeHexa::dNdr203(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr203[igauss][ishape][axis];
}
vvdouble& CShapeHexa::dNdr81(const uiint& igauss)
{
    return mvdNdr81[igauss];
}
vvdouble& CShapeHexa::dNdr82(const uiint& igauss)
{
    return mvdNdr82[igauss];
}
vvdouble& CShapeHexa::dNdr201(const uiint& igauss)
{
    return mvdNdr201[igauss];
}
vvdouble& CShapeHexa::dNdr202(const uiint& igauss)
{
    return mvdNdr202[igauss];
}
vvdouble& CShapeHexa::dNdr203(const uiint& igauss)
{
    return mvdNdr203[igauss];
}
double* CShapeHexa::Weight3d(const uiint& integNum1D)
{
    switch(integNum1D) {
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
double& CShapeHexa::Weight3dpt1()
{
    return mW3d1[0];
}
double& CShapeHexa::Weight3dpt2(const uiint& igauss)
{
    return mW3d2[igauss];
}
double& CShapeHexa::Weight3dpt3(const uiint& igauss)
{
    return mW3d3[igauss];
}
double& CShapeHexa::Weight3dpt4(const uiint& igauss)
{
    return mW3d4[igauss];
}
double& CShapeHexa::Weight3dpt5(const uiint& igauss)
{
    return mW3d5[igauss];
}
void CShapeHexa::Calc_dNdx8(const uiint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(8);
    uiint dof(3);
    switch(numOfInteg) {
    case(1):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr81, mpJacobi81, vNode, mvdNdx81, mv_detJ81);
        break;
    case(8):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr82, mpJacobi82, vNode, mvdNdx82, mv_detJ82);
        break;
    }
}
void CShapeHexa::Calc_dNdx20(const uiint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(20);
    uiint dof(3);
    switch(numOfInteg) {
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
double& CShapeHexa::detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, const uiint& ishape)
{
    switch(elemType) {
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
void CShapeHexa::setupIntegValue20()
{
    uiint ishape;
    for(ishape=0; ishape < 20; ishape++) mvIntegValue20[ishape]= 0.0;
    uiint igauss;
    for(igauss=0; igauss < 27; igauss++) {
        for(ishape=0; ishape < 20; ishape++) {
            mvIntegValue20[ishape] += mvN203[igauss][ishape];
        };
    };
    for(ishape=0; ishape < 20; ishape++) mvIntegValue20[ishape] /= 27.0;
}
void CShapeHexa::setupIntegValue8()
{
    uiint ishape;
    for(ishape=0; ishape < 8; ishape++) {
        mvIntegValue8[ishape]= 0.125;
    };
}
