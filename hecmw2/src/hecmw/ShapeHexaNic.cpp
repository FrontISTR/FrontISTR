/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeHexaNic.cpp
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
#include "Jacobian.h"
#include "ShapeHexaNic.h"
using namespace pmw;
CShapeHexaNic::CShapeHexaNic()
{
    mGzi1[0][0]= 0.0;
    mGzi1[0][1]= 0.0;
    mGzi1[0][2]= 0.0;
    mGzi8[0][0]= -0.577350269189626;
    mGzi8[0][1]= -0.577350269189626;
    mGzi8[0][2]= -0.577350269189626;
    mGzi8[1][0]=  0.577350269189626;
    mGzi8[1][1]= -0.577350269189626;
    mGzi8[1][2]= -0.577350269189626;
    mGzi8[2][0]= -0.577350269189626;
    mGzi8[2][1]=  0.577350269189626;
    mGzi8[2][2]= -0.577350269189626;
    mGzi8[3][0]=  0.577350269189626;
    mGzi8[3][1]=  0.577350269189626;
    mGzi8[3][2]= -0.577350269189626;
    mGzi8[4][0]= -0.577350269189626;
    mGzi8[4][1]= -0.577350269189626;
    mGzi8[4][2]=  0.577350269189626;
    mGzi8[5][0]=  0.577350269189626;
    mGzi8[5][1]= -0.577350269189626;
    mGzi8[5][2]=  0.577350269189626;
    mGzi8[6][0]= -0.577350269189626;
    mGzi8[6][1]=  0.577350269189626;
    mGzi8[6][2]=  0.577350269189626;
    mGzi8[7][0]=  0.577350269189626;
    mGzi8[7][1]=  0.577350269189626;
    mGzi8[7][2]=  0.577350269189626;
    mGzi27[0][0]=  -0.7745966692414830;
    mGzi27[0][1]=  -0.7745966692414830;
    mGzi27[0][2]=  -0.7745966692414830;
    mGzi27[1][0]=   0.0;
    mGzi27[1][1]=  -0.7745966692414830;
    mGzi27[1][2]=  -0.7745966692414830;
    mGzi27[2][0]=   0.7745966692414830;
    mGzi27[2][1]=  -0.7745966692414830;
    mGzi27[2][2]=  -0.7745966692414830;
    mGzi27[3][0]=  -0.7745966692414830;
    mGzi27[3][1]=   0.0;
    mGzi27[3][2]=  -0.7745966692414830;
    mGzi27[4][0]=   0.0;
    mGzi27[4][1]=   0.0;
    mGzi27[4][2]=  -0.7745966692414830;
    mGzi27[5][0]=   0.7745966692414830;
    mGzi27[5][1]=   0.0;
    mGzi27[5][2]=  -0.7745966692414830;
    mGzi27[6][0]=  -0.7745966692414830;
    mGzi27[6][1]=   0.7745966692414830;
    mGzi27[6][2]=  -0.7745966692414830;
    mGzi27[7][0]=   0.0;
    mGzi27[7][1]=   0.7745966692414830;
    mGzi27[7][2]=  -0.7745966692414830;
    mGzi27[8][0]=   0.7745966692414830;
    mGzi27[8][1]=   0.7745966692414830;
    mGzi27[8][2]=  -0.7745966692414830;
    mGzi27[9][0]=  -0.7745966692414830;
    mGzi27[9][1]=  -0.7745966692414830;
    mGzi27[9][2]=   0.0;
    mGzi27[10][0]=  0.0;
    mGzi27[10][1]= -0.7745966692414830;
    mGzi27[10][2]=  0.0;
    mGzi27[11][0]=  0.7745966692414830;
    mGzi27[11][1]= -0.7745966692414830;
    mGzi27[11][2]=  0.0;
    mGzi27[12][0]= -0.7745966692414830;
    mGzi27[12][1]=  0.0;
    mGzi27[12][2]=  0.0;
    mGzi27[13][0]=  0.0;
    mGzi27[13][1]=  0.0;
    mGzi27[13][2]=  0.0;
    mGzi27[14][0]=  0.7745966692414830;
    mGzi27[14][1]=  0.0;
    mGzi27[14][2]=  0.0;
    mGzi27[15][0]= -0.7745966692414830;
    mGzi27[15][1]=  0.7745966692414830;
    mGzi27[15][2]=  0.0;
    mGzi27[16][0]=  0.0;
    mGzi27[16][1]=  0.7745966692414830;
    mGzi27[16][2]=  0.0;
    mGzi27[17][0]=  0.7745966692414830;
    mGzi27[17][1]=  0.7745966692414830;
    mGzi27[17][2]=  0.0;
    mGzi27[18][0]= -0.7745966692414830;
    mGzi27[18][1]= -0.7745966692414830;
    mGzi27[18][2]=  0.7745966692414830;
    mGzi27[19][0]=  0.0;
    mGzi27[19][1]= -0.7745966692414830;
    mGzi27[19][2]=  0.7745966692414830;
    mGzi27[20][0]=  0.7745966692414830;
    mGzi27[20][1]= -0.7745966692414830;
    mGzi27[20][2]=  0.7745966692414830;
    mGzi27[21][0]= -0.7745966692414830;
    mGzi27[21][1]=  0.0;
    mGzi27[21][2]=  0.7745966692414830;
    mGzi27[22][0]=  0.0;
    mGzi27[22][1]=  0.0;
    mGzi27[22][2]=  0.7745966692414830;
    mGzi27[23][0]=  0.7745966692414830;
    mGzi27[23][1]=  0.0;
    mGzi27[23][2]=  0.7745966692414830;
    mGzi27[24][0]= -0.7745966692414830;
    mGzi27[24][1]=  0.7745966692414830;
    mGzi27[24][2]=  0.7745966692414830;
    mGzi27[25][0]=  0.0;
    mGzi27[25][1]=  0.7745966692414830;
    mGzi27[25][2]=  0.7745966692414830;
    mGzi27[26][0]=  0.7745966692414830;
    mGzi27[26][1]=  0.7745966692414830;
    mGzi27[26][2]=  0.7745966692414830;
    mW1[0]= 8.0;
    for(uiint i=0; i < 8; i++) mW8[i]= 1.0;
    mW27[0]=  0.171467764060357;
    mW27[1]=  0.274348422496571;
    mW27[2]=  0.171467764060357;
    mW27[3]=  0.274348422496571;
    mW27[4]=  0.438957475994513;
    mW27[5]=  0.274348422496571;
    mW27[6]=  0.171467764060357;
    mW27[7]=  0.274348422496571;
    mW27[8]=  0.171467764060357;
    mW27[9]=  0.274348422496571;
    mW27[10]= 0.438957475994513;
    mW27[11]= 0.274348422496571;
    mW27[12]= 0.438957475994513;
    mW27[13]= 0.702331961591221;
    mW27[14]= 0.438957475994513;
    mW27[15]= 0.274348422496571;
    mW27[16]= 0.438957475994513;
    mW27[17]= 0.274348422496571;
    mW27[18]= 0.171467764060357;
    mW27[19]= 0.274348422496571;
    mW27[20]= 0.171467764060357;
    mW27[21]= 0.274348422496571;
    mW27[22]= 0.438957475994513;
    mW27[23]= 0.274348422496571;
    mW27[24]= 0.171467764060357;
    mW27[25]= 0.274348422496571;
    mW27[26]= 0.171467764060357;
    mvIntegNum.resize(3);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 8;
    mvIntegNum[2]= 27;
    mvIntegNum1D.resize(3);
    mvIntegNum1D[0]= 1;
    mvIntegNum1D[1]= 2;
    mvIntegNum1D[2]= 3;
    uiint numOfIntg,numOfShape,dof;
    numOfIntg=1;
    numOfShape=11;
    dof=3;
    ResizeShape(mvN111, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr111, numOfIntg,numOfShape,dof);
    numOfIntg=8;
    numOfShape=11;
    dof=3;
    ResizeShape(mvN118, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr118, numOfIntg,numOfShape,dof);
    numOfIntg=27;
    numOfShape=11;
    dof=3;
    ResizeShape(mvN1127, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr1127, numOfIntg,numOfShape,dof);
    mpJacobi111= new CJacobian;
    mpJacobi118= new CJacobian;
    mpJacobi1127= new CJacobian;
    mpJacobi111->setupRegion(1,11);
    mpJacobi118->setupRegion(8,11);
    mpJacobi1127->setupRegion(27,11);
    numOfIntg=1;
    numOfShape=11;
    dof=3;
    ResizeDeriv(mvdNdx111, numOfIntg,numOfShape,dof);
    numOfIntg=8;
    numOfShape=11;
    dof=3;
    ResizeDeriv(mvdNdx118, numOfIntg,numOfShape,dof);
    numOfIntg=27;
    numOfShape=11;
    dof=3;
    ResizeDeriv(mvdNdx1127, numOfIntg,numOfShape,dof);
    numOfIntg= 1;
    numOfShape= 11;
    Resize_detJ(mv_detJ111, numOfIntg);
    numOfIntg= 8;
    numOfShape= 11;
    Resize_detJ(mv_detJ118, numOfIntg);
    numOfIntg= 27;
    numOfShape= 11;
    Resize_detJ(mv_detJ1127,numOfIntg);
    setupShapeFunction();
    setupShapeDeriv();
}
CShapeHexaNic::~CShapeHexaNic()
{
    ;
}
void CShapeHexaNic::setupShapeFunction()
{
    uiint igauss;
    double r,s,t;
    uiint numOfIntg;
    r= mGzi1[0][0];
    s= mGzi1[0][1];
    t= mGzi1[0][2];
    ShapeFunction11(mvN111, 0, r,s,t);
    numOfIntg= 8;
    for(igauss=0; igauss < numOfIntg; igauss++) {
        r= mGzi8[igauss][0];
        s= mGzi8[igauss][1];
        t= mGzi8[igauss][2];
        ShapeFunction11(mvN118, igauss, r,s,t);
    };
    numOfIntg= 27;
    for(igauss=0; igauss < numOfIntg; igauss++) {
        r= mGzi27[igauss][0];
        s= mGzi27[igauss][1];
        t= mGzi27[igauss][2];
        ShapeFunction11(mvN1127, igauss, r,s,t);
    }
}
void CShapeHexaNic::setupShapeDeriv()
{
    uiint igauss;
    double r,s,t;
    uiint numOfIntg;
    r= mGzi1[0][0];
    s= mGzi1[0][1];
    t= mGzi1[0][2];
    ShapeDeriv11(mvdNdr111, 0, r,s,t);
    numOfIntg= 8;
    for(igauss=0; igauss < numOfIntg; igauss++) {
        r= mGzi8[igauss][0];
        s= mGzi8[igauss][1];
        t= mGzi8[igauss][2];
        ShapeDeriv11(mvdNdr118, igauss, r,s,t);
    };
    numOfIntg= 27;
    for(igauss=0; igauss< numOfIntg; igauss++) {
        r= mGzi27[igauss][0];
        s= mGzi27[igauss][1];
        t= mGzi27[igauss][2];
        ShapeDeriv11(mvdNdr1127, igauss, r,s,t);
    }
}
void CShapeHexaNic::ShapeFunction11(vvdouble& N, const uiint& igauss, const double& r, const double& s, const double& t)
{
    N[igauss][0]=  0.125*(1.0-r)*(1.0-s)*(1.0-t);
    N[igauss][1]=  0.125*(1.0+r)*(1.0-s)*(1.0-t);
    N[igauss][2]=  0.125*(1.0+r)*(1.0+s)*(1.0-t);
    N[igauss][3]=  0.125*(1.0-r)*(1.0+s)*(1.0-t);
    N[igauss][4]=  0.125*(1.0-r)*(1.0-s)*(1.0+t);
    N[igauss][5]=  0.125*(1.0+r)*(1.0-s)*(1.0+t);
    N[igauss][6]=  0.125*(1.0+r)*(1.0+s)*(1.0+t);
    N[igauss][7]=  0.125*(1.0-r)*(1.0+s)*(1.0+t);
    N[igauss][8]=  1.0-r*r;
    N[igauss][9]=  1.0-s*s;
    N[igauss][10]= 1.0-t*t;
}
void CShapeHexaNic::ShapeDeriv11(vvvdouble& dNdr, const uiint& igauss, const double& r, const double& s, const double& t)
{
    dNdr[igauss][0][0]= -0.125*(1.0-s)*(1.0-t);
    dNdr[igauss][1][0]=  0.125*(1.0-s)*(1.0-t);
    dNdr[igauss][2][0]=  0.125*(1.0+s)*(1.0-t);
    dNdr[igauss][3][0]= -0.125*(1.0+s)*(1.0-t);
    dNdr[igauss][4][0]= -0.125*(1.0-s)*(1.0+t);
    dNdr[igauss][5][0]=  0.125*(1.0-s)*(1.0+t);
    dNdr[igauss][6][0]=  0.125*(1.0+s)*(1.0+t);
    dNdr[igauss][7][0]= -0.125*(1.0+s)*(1.0+t);
    dNdr[igauss][8][0]= -2.0*r;
    dNdr[igauss][9][0]=  0.0;
    dNdr[igauss][10][0]= 0.0;
    dNdr[igauss][0][1]= -0.125*(1.0-r)*(1.0-t);
    dNdr[igauss][1][1]= -0.125*(1.0+r)*(1.0-t);
    dNdr[igauss][2][1]=  0.125*(1.0+r)*(1.0-t);
    dNdr[igauss][3][1]=  0.125*(1.0-r)*(1.0-t);
    dNdr[igauss][4][1]= -0.125*(1.0-r)*(1.0+t);
    dNdr[igauss][5][1]= -0.125*(1.0+r)*(1.0+t);
    dNdr[igauss][6][1]=  0.125*(1.0+r)*(1.0+t);
    dNdr[igauss][7][1]=  0.125*(1.0-r)*(1.0+t);
    dNdr[igauss][8][1]=  0.0;
    dNdr[igauss][9][1]= -2.0*s;
    dNdr[igauss][10][1]= 0.0;
    dNdr[igauss][0][2]= -0.125*(1.0-r)*(1.0-s);
    dNdr[igauss][1][2]= -0.125*(1.0+r)*(1.0-s);
    dNdr[igauss][2][2]= -0.125*(1.0+r)*(1.0+s);
    dNdr[igauss][3][2]= -0.125*(1.0-r)*(1.0+s);
    dNdr[igauss][4][2]=  0.125*(1.0-r)*(1.0-s);
    dNdr[igauss][5][2]=  0.125*(1.0+r)*(1.0-s);
    dNdr[igauss][6][2]=  0.125*(1.0+r)*(1.0+s);
    dNdr[igauss][7][2]=  0.125*(1.0-r)*(1.0+s);
    dNdr[igauss][8][2]=  0.0;
    dNdr[igauss][9][2]=  0.0;
    dNdr[igauss][10][2]= -2.0*t;
}
double& CShapeHexaNic::N111(const uiint& igauss, const uiint& ishape)
{
    return mvN111[igauss][ishape];
}
double& CShapeHexaNic::N118(const uiint& igauss, const uiint& ishape)
{
    return mvN118[igauss][ishape];
}
double& CShapeHexaNic::N1127(const uiint& igauss, const uiint& ishape)
{
    return mvN1127[igauss][ishape];
}
vdouble& CShapeHexaNic::N111(const uiint& igauss)
{
    return mvN111[igauss];
}
vdouble& CShapeHexaNic::N118(const uiint& igauss)
{
    return mvN118[igauss];
}
vdouble& CShapeHexaNic::N1127(const uiint& igauss)
{
    return mvN1127[igauss];
}
double& CShapeHexaNic::dNdr111(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr111[igauss][ishape][axis];
}
double& CShapeHexaNic::dNdr118(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr118[igauss][ishape][axis];
}
double& CShapeHexaNic::dNdr1127(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr1127[igauss][ishape][axis];
}
vvdouble& CShapeHexaNic::dNdr111(const uiint& igauss)
{
    return mvdNdr111[igauss];
}
vvdouble& CShapeHexaNic::dNdr118(const uiint& igauss)
{
    return mvdNdr118[igauss];
}
vvdouble& CShapeHexaNic::dNdr1127(const uiint& igauss)
{
    return mvdNdr1127[igauss];
}
double* CShapeHexaNic::Weight(const uiint& integNum)
{
    switch(integNum) {
    case(1):
        return mW1;
        break;
    case(8):
        return mW8;
        break;
    case(27):
        return mW27;
        break;
    default:
        break;
    }
}
double& CShapeHexaNic::Weight_pt1()
{
    return mW1[0];
}
double& CShapeHexaNic::Weight_pt8(const uiint& igauss)
{
    return mW8[igauss];
}
double& CShapeHexaNic::Weight_pt27(const uiint& igauss)
{
    return mW27[igauss];
}
void CShapeHexaNic::Calc_dNdx11(const uiint& numOfInteg, CElement* pElement)
{
}
double& CShapeHexaNic::detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss)
{
    switch(elemType) {
    case(ElementType::HexaNic):
        if(numOfInteg==1)
            return mv_detJ111[igauss];
        if(numOfInteg==8)
            return mv_detJ118[igauss];
        if(numOfInteg==27)
            return mv_detJ1127[igauss];
        break;
    }
}
