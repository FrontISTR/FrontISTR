/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeTetra.h
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
#include "TypeDef.h"
#include "ShapeFunctionBase.h"
namespace pmw
{
#ifndef _SHAPETETRA_H
#define	_SHAPETETRA_H
class CShapeTetra:public CShapeFunctionBase
{
private:
    CShapeTetra();
public:
    static CShapeTetra* Instance() {
        static CShapeTetra moShapeTetra;
        return &moShapeTetra;
    }
    virtual ~CShapeTetra();
private:
    double mGzi1[1][3];
    double mGzi4[4][3];
    double mGzi15[15][3];
    double mW1[1];
    double mW4[4];
    double mW15[15];
    vvdouble mvN41;
    vvdouble mvN101;
    vvdouble mvN104;
    vvdouble mvN1015;
    vvvdouble mvdNdr41;
    vvvdouble mvdNdr101;
    vvvdouble mvdNdr104;
    vvvdouble mvdNdr1015;
    CJacobian *mpJacobi41;
    CJacobian *mpJacobi101;
    CJacobian *mpJacobi104;
    CJacobian *mpJacobi1015;
    vvvdouble mvdNdx41;
    vvvdouble mvdNdx101;
    vvvdouble mvdNdx104;
    vvvdouble mvdNdx1015;
    vdouble mv_detJ41;
    vdouble mv_detJ101;
    vdouble mv_detJ104;
    vdouble mv_detJ1015;
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void ShapeFunction4(vvdouble& N, const uiint& igauss,
                        const double& L1, const double& L2, const double& L3);
    void ShapeFunction10(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& L3);
    void ShapeDeriv4(vvvdouble& dNdr, const uiint& igauss);
    void ShapeDeriv10(vvvdouble& dNdr, const uiint& igauss,
                      const double& L1, const double& L2, const double& L3);
    vdouble mvIntegValue10;
    vdouble mvIntegValue4;
    void setupIntegValue10();
    void setupIntegValue4();
public:
    double& N41(const uiint& igauss, const uiint& ishape);
    double& N101(const uiint& igauss, const uiint& ishape);
    double& N104(const uiint& igauss, const uiint& ishape);
    double& N1015(const uiint& igauss, const uiint& ishape);
    vdouble& N41(const uiint& igauss) {
        return mvN41[igauss];
    }
    vdouble& N101(const uiint& igauss) {
        return mvN101[igauss];
    }
    vdouble& N104(const uiint& igauss) {
        return mvN104[igauss];
    }
    vdouble& N1015(const uiint& igauss) {
        return mvN1015[igauss];
    }
    vvdouble& N41() {
        return mvN41;
    }
    vvdouble& N101() {
        return mvN101;
    }
    vvdouble& N104() {
        return mvN104;
    }
    vvdouble& N1015() {
        return mvN1015;
    }
    double& dNdr41(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr101(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr104(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr1015(const uiint& igauss, const uiint& ishape, const uiint& axis);
    vvdouble& dNdr41(const uiint& igauss) {
        return mvdNdr41[igauss];
    }
    vvdouble& dNdr101(const uiint& igauss) {
        return mvdNdr101[igauss];
    }
    vvdouble& dNdr104(const uiint& igauss) {
        return mvdNdr104[igauss];
    }
    vvdouble& dNdr1015(const uiint& igauss) {
        return mvdNdr1015[igauss];
    }
    vvvdouble& dNdr41() {
        return mvdNdr41;
    }
    vvvdouble& dNdr101() {
        return mvdNdr101;
    }
    vvvdouble& dNdr104() {
        return mvdNdr104;
    }
    vvvdouble& dNdr1015() {
        return mvdNdr1015;
    }
    double* Weight(const uiint& integNum);
    double& Weight_pt1();
    double& Weight_pt4(const uiint& igauss);
    double& Weight_pt15(const uiint& igauss);
    void Calc_dNdx4(const uiint& numOfInteg, CElement *pElement);
    void Calc_dNdx10(const uiint& numOfInteg, CElement *pElement);
    double& dNdx41(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return  mvdNdx41[igauss][ishape][axis];
    }
    double& dNdx101(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx101[igauss][ishape][axis];
    }
    double& dNdx104(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx104[igauss][ishape][axis];
    }
    double& dNdx1015(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx1015[igauss][ishape][axis];
    }
    vvdouble& dNdx41(const uiint& igauss) {
        return  mvdNdx41[igauss];
    }
    vvdouble& dNdx101(const uiint& igauss) {
        return mvdNdx101[igauss];
    }
    vvdouble& dNdx104(const uiint& igauss) {
        return mvdNdx104[igauss];
    }
    vvdouble& dNdx1015(const uiint& igauss) {
        return mvdNdx1015[igauss];
    }
    vvvdouble& dNdx41() {
        return   mvdNdx41;
    }
    vvvdouble& dNdx101() {
        return  mvdNdx101;
    }
    vvvdouble& dNdx104() {
        return  mvdNdx104;
    }
    vvvdouble& dNdx1015() {
        return mvdNdx1015;
    }
    double& detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss);
    double& detJ41(const uiint& igauss) {
        return mv_detJ41[igauss];
    }
    double& detJ101(const uiint& igauss) {
        return mv_detJ101[igauss];
    }
    double& detJ104(const uiint& igauss) {
        return mv_detJ104[igauss];
    }
    double& detJ1015(const uiint& igauss) {
        return mv_detJ1015[igauss];
    }
    vdouble& detJ41() {
        return   mv_detJ41;
    }
    vdouble& detJ101() {
        return  mv_detJ101;
    }
    vdouble& detJ104() {
        return  mv_detJ104;
    }
    vdouble& detJ1015() {
        return mv_detJ1015;
    }
    vdouble& getIntegValue10() {
        return mvIntegValue10;
    }
    vdouble& getIntegValue4() {
        return mvIntegValue4;
    }
    double& getIntegValue10(const uiint& ishape) {
        return mvIntegValue10[ishape];
    }
    double& getIntegValue4(const uiint& ishape) {
        return mvIntegValue4[ishape];
    }
};
#endif	/* _SHAPETETRA_H */
}
