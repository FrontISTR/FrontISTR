/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapePrism.h
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
#ifndef _SHAPEPRISM_H
#define	_SHAPEPRISM_H
class CShapePrism:public CShapeFunctionBase
{
private:
    CShapePrism();
public:
    static CShapePrism* Instance() {
        static CShapePrism moShapePrism;
        return &moShapePrism;
    }
    virtual ~CShapePrism();
private:
    double mGzi2[2][3];
    double mGzi6[6][3];
    double mGzi9[9][3];
    double mGzi18[18][3];
    double mW2[2];
    double mW6[6];
    double mW9[9];
    double mW18[18];
    vvdouble mvN62;
    vvdouble mvN156;
    vvdouble mvN159;
    vvdouble mvN1518;
    vvvdouble mvdNdr62;
    vvvdouble mvdNdr156;
    vvvdouble mvdNdr159;
    vvvdouble mvdNdr1518;
    CJacobian *mpJacobi62;
    CJacobian *mpJacobi156;
    CJacobian *mpJacobi159;
    CJacobian *mpJacobi1518;
    vvvdouble mvdNdx62;
    vvvdouble mvdNdx156;
    vvvdouble mvdNdx159;
    vvvdouble mvdNdx1518;
    vdouble mv_detJ62;
    vdouble mv_detJ156;
    vdouble mv_detJ159;
    vdouble mv_detJ1518;
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void ShapeFunction6(vvdouble& N, const uiint& igauss,
                        const double& L1, const double& L2, const double& gzi);
    void ShapeFunction15(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& gzi);
    void ShapeDeriv6(vvvdouble& dNdr, const uiint& igauss,
                     const double& L1, const double& L2, const double& gzi);
    void ShapeDeriv15(vvvdouble& dNdr, const uiint& igauss,
                      const double& L1, const double& L2, const double& gzi);
    vdouble mvIntegValue15;
    vdouble mvIntegValue6;
    void setupIntegValue15();
    void setupIntegValue6();
public:
    double& N62(const uiint& igauss, const uiint& ishape);
    double& N156(const uiint& igauss, const uiint& ishape);
    double& N159(const uiint& igauss, const uiint& ishape);
    double& N1518(const uiint& igauss, const uiint& ishape);
    vdouble& N62(const uiint& igauss) {
        return mvN62[igauss];
    }
    vdouble& N156(const uiint& igauss) {
        return mvN156[igauss];
    }
    vdouble& N159(const uiint& igauss) {
        return mvN159[igauss];
    }
    vdouble& N1518(const uiint& igauss) {
        return mvN1518[igauss];
    }
    vvdouble& N62() {
        return mvN62;
    }
    vvdouble& N156() {
        return mvN156;
    }
    vvdouble& N159() {
        return mvN159;
    }
    vvdouble& N1518() {
        return mvN1518;
    }
    double& dNdr62(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr156(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr159(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr1518(const uiint& igauss, const uiint& ishape, const uiint& axis);
    vvdouble& dNdr62(const uiint& igauss) {
        return mvdNdr62[igauss];
    }
    vvdouble& dNdr156(const uiint& igauss) {
        return mvdNdr156[igauss];
    }
    vvdouble& dNdr159(const uiint& igauss) {
        return mvdNdr159[igauss];
    }
    vvdouble& dNdr1518(const uiint& igauss) {
        return mvdNdr1518[igauss];
    }
    vvvdouble& dNdr62() {
        return mvdNdr62;
    }
    vvvdouble& dNdr156() {
        return mvdNdr156;
    }
    vvvdouble& dNdr159() {
        return mvdNdr159;
    }
    vvvdouble& dNdr1518() {
        return mvdNdr1518;
    }
    double* Weight(const uiint& integNum);
    double& Weight_pt2(const uiint& igauss);
    double& Weight_pt6(const uiint& igauss);
    double& Weight_pt9(const uiint& igauss);
    double& Weight_pt18(const uiint& igauss);
    void Calc_dNdx6(const uiint& numOfInteg, CElement *pElement);
    void Calc_dNdx15(const uiint& numOfInteg, CElement *pElement);
    double& dNdx62(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return  mvdNdx62[igauss][ishape][axis];
    }
    double& dNdx156(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx156[igauss][ishape][axis];
    }
    double& dNdx159(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx159[igauss][ishape][axis];
    }
    double& dNdx1518(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx1518[igauss][ishape][axis];
    }
    vvdouble& dNdx62(const uiint& igauss) {
        return  mvdNdx62[igauss];
    }
    vvdouble& dNdx156(const uiint& igauss) {
        return mvdNdx156[igauss];
    }
    vvdouble& dNdx159(const uiint& igauss) {
        return mvdNdx159[igauss];
    }
    vvdouble& dNdx1518(const uiint& igauss) {
        return mvdNdx1518[igauss];
    }
    vvvdouble& dNdx62() {
        return   mvdNdx62;
    }
    vvvdouble& dNdx156() {
        return  mvdNdx156;
    }
    vvvdouble& dNdx159() {
        return  mvdNdx159;
    }
    vvvdouble& dNdx1518() {
        return mvdNdx1518;
    }
    double& detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss);
    double& detJ62(const uiint& igauss) {
        return mv_detJ62[igauss];
    }
    double& detJ156(const uiint& igauss) {
        return mv_detJ156[igauss];
    }
    double& detJ159(const uiint& igauss) {
        return mv_detJ159[igauss];
    }
    double& detJ1518(const uiint& igauss) {
        return mv_detJ1518[igauss];
    }
    vdouble& detJ62() {
        return   mv_detJ62;
    }
    vdouble& detJ156() {
        return  mv_detJ156;
    }
    vdouble& detJ159() {
        return  mv_detJ159;
    }
    vdouble& detJ1518() {
        return mv_detJ1518;
    }
    vdouble& getIntegValue15() {
        return mvIntegValue15;
    }
    vdouble& getIntegValue6() {
        return mvIntegValue6;
    }
    double& getIntegValue15(const uiint& ishape) {
        return mvIntegValue15[ishape];
    }
    double& getIntegValue6(const uiint& ishape) {
        return mvIntegValue6[ishape];
    }
};
#endif	/* _SHAPEPRISM_H */
}
