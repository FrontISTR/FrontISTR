/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeTetra.cpp
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
#include "ShapeTetra.h"
using namespace pmw;
CShapeTetra::CShapeTetra()
{
    mGzi1[0][0]= 0.25;
    mGzi1[0][1]= 0.25;
    mGzi1[0][2]= 0.25;
    mGzi4[0][0]= 0.138196601125011;
    mGzi4[0][1]= 0.138196601125011;
    mGzi4[0][2]= 0.138196601125011;
    mGzi4[1][0]= 0.585410196624968;
    mGzi4[1][1]= 0.138196601125011;
    mGzi4[1][2]= 0.138196601125011;
    mGzi4[2][0]= 0.138196601125011;
    mGzi4[2][1]= 0.585410196624968;
    mGzi4[2][2]= 0.138196601125011;
    mGzi4[3][0]= 0.138196601125011;
    mGzi4[3][1]= 0.138196601125011;
    mGzi4[3][2]= 0.585410196624968;
    mGzi15[0][0]=  0.25;
    mGzi15[0][1]=  0.25;
    mGzi15[0][2]=  0.25;
    mGzi15[1][0]=  0.091971078052723;
    mGzi15[1][1]=  0.091971078052723;
    mGzi15[1][2]=  0.091971078052723;
    mGzi15[2][0]=  0.724086765841831;
    mGzi15[2][1]=  0.091971078052723;
    mGzi15[2][2]=  0.091971078052723;
    mGzi15[3][0]=  0.091971078052723;
    mGzi15[3][1]=  0.724086765841831;
    mGzi15[3][2]=  0.091971078052723;
    mGzi15[4][0]=  0.091971078052723;
    mGzi15[4][1]=  0.091971078052723;
    mGzi15[4][2]=  0.724086765841831;
    mGzi15[5][0]=  0.319793627829630;
    mGzi15[5][1]=  0.319793627829630;
    mGzi15[5][2]=  0.319793627829630;
    mGzi15[6][0]=  0.040619116511110;
    mGzi15[6][1]=  0.319793627829630;
    mGzi15[6][2]=  0.319793627829630;
    mGzi15[7][0]=  0.319793627829630;
    mGzi15[7][1]=  0.040619116511110;
    mGzi15[7][2]=  0.319793627829630;
    mGzi15[8][0]=  0.319793627829630;
    mGzi15[8][1]=  0.319793627829630;
    mGzi15[8][2]=  0.040619116511110;
    mGzi15[9][0]=  0.056350832689629;
    mGzi15[9][1]=  0.056350832689629;
    mGzi15[9][2]=  0.443649167310371;
    mGzi15[10][0]= 0.443649167310371;
    mGzi15[10][1]= 0.056350832689629;
    mGzi15[10][2]= 0.056350832689629;
    mGzi15[11][0]= 0.443649167310371;
    mGzi15[11][1]= 0.443649167310371;
    mGzi15[11][2]= 0.056350832689629;
    mGzi15[12][0]= 0.056350832689629;
    mGzi15[12][1]= 0.443649167310371;
    mGzi15[12][2]= 0.443649167310371;
    mGzi15[13][0]= 0.056350832689629;
    mGzi15[13][1]= 0.443649167310371;
    mGzi15[13][2]= 0.056350832689629;
    mGzi15[14][0]= 0.443649167310371;
    mGzi15[14][1]= 0.056350832689629;
    mGzi15[14][2]= 0.443649167310371;
    mW1[0]= 0.166666666666667;
    mW4[0]= 0.041666666666667;
    mW4[1]= 0.041666666666667;
    mW4[2]= 0.041666666666667;
    mW4[3]= 0.041666666666667;
    mW15[0]=  0.019753086419753;
    mW15[1]=  0.011989513963170;
    mW15[2]=  0.011989513963170;
    mW15[3]=  0.011989513963170;
    mW15[4]=  0.011989513963170;
    mW15[5]=  0.011511367871045;
    mW15[6]=  0.011511367871045;
    mW15[7]=  0.011511367871045;
    mW15[8]=  0.011511367871045;
    mW15[9]=  0.008818342151675;
    mW15[10]= 0.008818342151675;
    mW15[11]= 0.008818342151675;
    mW15[12]= 0.008818342151675;
    mW15[13]= 0.008818342151675;
    mW15[14]= 0.008818342151675;
    mvIntegNum.resize(3);
    mvIntegNum[0]= 1;
    mvIntegNum[1]= 4;
    mvIntegNum[2]= 15;
    uiint numOfIntg,numOfShape,dof;
    numOfIntg=1;
    numOfShape=4;
    dof=3;
    ResizeShape(mvN41, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr41, numOfIntg,numOfShape,dof);
    numOfIntg=1;
    numOfShape=10;
    dof=3;
    ResizeShape(mvN101, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr101, numOfIntg,numOfShape,dof);
    numOfIntg=4;
    numOfShape=10;
    dof=3;
    ResizeShape(mvN104, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr104, numOfIntg,numOfShape,dof);
    numOfIntg=15;
    numOfShape=10;
    dof=3;
    ResizeShape(mvN1015, numOfIntg,numOfShape);
    ResizeDeriv(mvdNdr1015, numOfIntg,numOfShape,dof);
    mpJacobi41= new CJacobian;
    mpJacobi101= new CJacobian;
    mpJacobi104= new CJacobian;
    mpJacobi1015= new CJacobian;
    mpJacobi41->setupRegion(1,4);
    mpJacobi101->setupRegion(1,10);
    mpJacobi104->setupRegion(4,10);
    mpJacobi1015->setupRegion(15,10);
    numOfIntg= 1;
    numOfShape= 4;
    dof= 3;
    ResizeDeriv(mvdNdx41, numOfIntg, numOfShape, dof);
    numOfIntg= 1;
    numOfShape= 10;
    dof= 3;
    ResizeDeriv(mvdNdx101, numOfIntg, numOfShape, dof);
    numOfIntg= 4;
    numOfShape= 10;
    dof= 3;
    ResizeDeriv(mvdNdx104, numOfIntg, numOfShape, dof);
    numOfIntg= 15;
    numOfShape= 10;
    dof= 3;
    ResizeDeriv(mvdNdx1015, numOfIntg, numOfShape, dof);
    numOfIntg= 1;
    numOfShape= 4;
    Resize_detJ(mv_detJ41, numOfIntg);
    numOfIntg= 1;
    numOfShape= 10;
    Resize_detJ(mv_detJ101, numOfIntg);
    numOfIntg= 4;
    numOfShape= 10;
    Resize_detJ(mv_detJ104, numOfIntg);
    numOfIntg= 15;
    numOfShape= 10;
    Resize_detJ(mv_detJ1015, numOfIntg);
    setupShapeFunction(mvN41, 1,4,mGzi1);
    setupShapeFunction(mvN101, 1,10,mGzi1);
    setupShapeFunction(mvN104, 4,10,mGzi4);
    setupShapeFunction(mvN1015, 15,10,mGzi15);
    uiint igauss= 0;
    {
        ShapeDeriv4(mvdNdr41,igauss);
    }
    setupShapeDeriv(mvdNdr101, 1,10,mGzi1);
    setupShapeDeriv(mvdNdr104, 4,10,mGzi4);
    setupShapeDeriv(mvdNdr1015, 15,10,mGzi15);
    mvIntegValue10.resize(10);
    mvIntegValue4.resize(4);
    setupIntegValue10();
    setupIntegValue4();
}
CShapeTetra::~CShapeTetra()
{
    ;
}
void CShapeTetra::setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3])
{
    double L1,L2,L3;
    uiint igauss;
    switch(numOfShape) {
    case(4):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            L1= Gzi[igauss][0];
            L2= Gzi[igauss][1];
            L3= Gzi[igauss][2];
            ShapeFunction4(N,igauss, L1,L2,L3);
        };
        break;
    case(10):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            L1= Gzi[igauss][0];
            L2= Gzi[igauss][1];
            L3= Gzi[igauss][2];
            ShapeFunction10(N,igauss, L1,L2,L3);
        };
        break;
    default:
        break;
    }
}
void CShapeTetra::setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3])
{
    double L1,L2,L3;
    uiint igauss;
    switch(numOfShape) {
    case(4):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            L1= Gzi[igauss][0];
            L2= Gzi[igauss][1];
            L3= Gzi[igauss][2];
            ShapeDeriv4(dNdr,igauss);
        };
        break;
    case(10):
        for(igauss=0; igauss< numOfIntg; igauss++) {
            L1= Gzi[igauss][0];
            L2= Gzi[igauss][1];
            L3= Gzi[igauss][2];
            ShapeDeriv10(dNdr,igauss, L1,L2,L3);
        };
        break;
    default:
        break;
    }
}
void CShapeTetra::ShapeFunction4(vvdouble& N, const uiint& igauss,
                                 const double& L1, const double& L2, const double& L3)
{
    N[igauss][0]= 1.0 - L1 - L2 - L3;
    N[igauss][1]= L1;
    N[igauss][2]= L2;
    N[igauss][3]= L3;
}
void CShapeTetra::ShapeFunction10(vvdouble& N, const uiint& igauss,
                                  const double& L1, const double& L2, const double& L3)
{
    double a= 1.0 - L1 - L2 - L3;
    N[igauss][0]= (2.0*a - 1.0)*a;
    N[igauss][1]= L1*(2.0*L1 - 1.0);
    N[igauss][2]= L2*(2.0*L2 - 1.0);
    N[igauss][3]= L3*(2.0*L3 - 1.0);
    N[igauss][4]= 4.0*L1*a;
    N[igauss][5]= 4.0*L1*L2;
    N[igauss][6]= 4.0*L2*a;
    N[igauss][7]= 4.0*L3*a;
    N[igauss][8]= 4.0*L1*L3;
    N[igauss][9]= 4.0*L2*L3;
}
void CShapeTetra::ShapeDeriv4(vvvdouble& dNdr, const uiint& igauss)
{
    dNdr[igauss][0][0] = -1.0;
    dNdr[igauss][1][0] =  1.0;
    dNdr[igauss][2][0] =  0.0;
    dNdr[igauss][3][0] =  0.0;
    dNdr[igauss][0][1] = -1.0;
    dNdr[igauss][1][1] =  0.0;
    dNdr[igauss][2][1] =  1.0;
    dNdr[igauss][3][1] =  0.0;
    dNdr[igauss][0][2] = -1.0;
    dNdr[igauss][1][2] =  0.0;
    dNdr[igauss][2][2] =  0.0;
    dNdr[igauss][3][2] =  1.0;
}
void CShapeTetra::ShapeDeriv10(vvvdouble& dNdr, const uiint& igauss,
                               const double& L1, const double& L2, const double& L3)
{
    double a= 1.0 - L1 - L2 - L3;
    dNdr[igauss][0][0]=  1.0 - 4.0*a;
    dNdr[igauss][1][0]=  4.0*L1 - 1.0;
    dNdr[igauss][2][0]=  0.0;
    dNdr[igauss][3][0]=  0.0;
    dNdr[igauss][4][0]=  4.0*(1.0 - 2.0*L1 - L2 - L3);
    dNdr[igauss][5][0]=  4.0*L2;
    dNdr[igauss][6][0]= -4.0*L2;
    dNdr[igauss][7][0]= -4.0*L3;
    dNdr[igauss][8][0]=  4.0*L3;
    dNdr[igauss][9][0]=  0.0;
    dNdr[igauss][0][1]=  1.0 - 4.0*a;
    dNdr[igauss][1][1]=  0.0;
    dNdr[igauss][2][1]=  4.0*L2-1.0;
    dNdr[igauss][3][1]=  0.0;
    dNdr[igauss][4][1]= -4.0*L1;
    dNdr[igauss][5][1]=  4.0*L1;
    dNdr[igauss][6][1]=  4.0*(1.0 - L1 - 2.0*L2 - L3);
    dNdr[igauss][7][1]= -4.0*L3;
    dNdr[igauss][8][1]=  0.0;
    dNdr[igauss][9][1]=  4.0*L3;
    dNdr[igauss][0][2]=  1.0 - 4.0*a;
    dNdr[igauss][1][2]=  0.0;
    dNdr[igauss][2][2]=  0.0;
    dNdr[igauss][3][2]=  4.0*L3 - 1.0;
    dNdr[igauss][4][2]= -4.0*L1;
    dNdr[igauss][5][2]=  0.0;
    dNdr[igauss][6][2]= -4.0*L2;
    dNdr[igauss][7][2]=  4.0*(1.0- L1 - L2 - 2.0*L3);
    dNdr[igauss][8][2]=  4.0*L1;
    dNdr[igauss][9][2]=  4.0*L2;
}
double& CShapeTetra::N41(const uiint& igauss, const uiint& ishape)
{
    return mvN41[igauss][ishape];
}
double& CShapeTetra::N101(const uiint& igauss, const uiint& ishape)
{
    return mvN101[igauss][ishape];
}
double& CShapeTetra::N104(const uiint& igauss, const uiint& ishape)
{
    return mvN104[igauss][ishape];
}
double& CShapeTetra::N1015(const uiint& igauss, const uiint& ishape)
{
    return mvN1015[igauss][ishape];
}
double& CShapeTetra::dNdr41(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr41[igauss][ishape][axis];
}
double& CShapeTetra::dNdr101(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr101[igauss][ishape][axis];
}
double& CShapeTetra::dNdr104(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr104[igauss][ishape][axis];
}
double& CShapeTetra::dNdr1015(const uiint& igauss, const uiint& ishape, const uiint& axis)
{
    return mvdNdr1015[igauss][ishape][axis];
}
double* CShapeTetra::Weight(const uiint& integNum)
{
    switch(integNum) {
    case(1):
        return mW1;
        break;
    case(4):
        return mW4;
        break;
    case(15):
        return mW15;
        break;
    default:
        break;
    }
}
double& CShapeTetra::Weight_pt1()
{
    return mW1[0];
}
double& CShapeTetra::Weight_pt4(const uiint& igauss)
{
    return mW4[igauss];
}
double& CShapeTetra::Weight_pt15(const uiint& igauss)
{
    return mW15[igauss];
}
void CShapeTetra::Calc_dNdx4(const uiint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(4);
    uiint dof(3);
    switch(numOfInteg) {
    case(1):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr41, mpJacobi41, vNode, mvdNdx41, mv_detJ41);
        break;
    default:
        break;
    }
}
void CShapeTetra::Calc_dNdx10(const uiint& numOfInteg, CElement *pElement)
{
    vector<CNode*> vNode= pElement->getNode();
    uiint numOfShape(10);
    uiint dof(3);
    switch(numOfInteg) {
    case(1):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr101, mpJacobi101, vNode, mvdNdx101, mv_detJ101);
        break;
    case(4):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr104, mpJacobi104, vNode, mvdNdx104, mv_detJ104);
        break;
    case(15):
        ShapeDerivXYZ(numOfInteg, numOfShape, dof, mvdNdr1015, mpJacobi1015, vNode, mvdNdx1015, mv_detJ1015);
        break;
    }
}
double& CShapeTetra::detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss)
{
    switch(elemType) {
    case(ElementType::Tetra):
        if(numOfInteg==1)
            return mv_detJ41[igauss];
        break;
    case(ElementType::Tetra2):
        if(numOfInteg==1)
            return mv_detJ101[igauss];
        if(numOfInteg==4)
            return mv_detJ104[igauss];
        if(numOfInteg==15)
            return mv_detJ1015[igauss];
        break;
    }
}
void CShapeTetra::setupIntegValue10()
{
    uiint ishape;
    for(ishape=0; ishape < 10; ishape++) mvIntegValue10[ishape]= 0.0;
    uiint igauss;
    for(igauss=0; igauss < 4; igauss++) {
        for(ishape=0; ishape < 10; ishape++) {
            mvIntegValue10[ishape] += mvN104[igauss][ishape];
        };
    };
    for(ishape=0; ishape < 4; ishape++) mvIntegValue10[ishape] /= 4.0;
}
void CShapeTetra::setupIntegValue4()
{
    uiint ishape;
    for(ishape=0; ishape < 4; ishape++) mvIntegValue4[ishape]= 0.25;
}
