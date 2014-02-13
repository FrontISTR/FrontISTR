/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeHexa.h
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
#ifndef _SHAPEHEXA_H
#define	_SHAPEHEXA_H
class CShapeHexa:public CShapeFunctionBase
{
private:
    CShapeHexa();
public:
    virtual ~CShapeHexa();
    static CShapeHexa* Instance() {
        static CShapeHexa moShapeHexa;
        return &moShapeHexa;
    }
private:
    double mGzi1[1];
    double mGzi2[2];
    double mGzi3[3];
    double mGzi4[4];
    double mGzi5[5];
    double mW1[1];
    double mW2[2];
    double mW3[3];
    double mW4[4];
    double mW5[5];
    double mW3d1[1];
    double mW3d2[8];
    double mW3d3[27];
    double mW3d4[64];
    double mW3d5[125];
    vvdouble mvN81;
    vvdouble mvN82;
    vvdouble mvN201;
    vvdouble mvN202;
    vvdouble mvN203;
    vvvdouble mvdNdr81;
    vvvdouble mvdNdr82;
    vvvdouble mvdNdr201;
    vvvdouble mvdNdr202;
    vvvdouble mvdNdr203;
    CJacobian *mpJacobi81;
    CJacobian *mpJacobi82;
    CJacobian *mpJacobi201;
    CJacobian *mpJacobi202;
    CJacobian *mpJacobi203;
    vvvdouble mvdNdx81;
    vvvdouble mvdNdx82;
    vvvdouble mvdNdx201;
    vvvdouble mvdNdx202;
    vvvdouble mvdNdx203;
    vdouble mv_detJ81;
    vdouble mv_detJ82;
    vdouble mv_detJ201;
    vdouble mv_detJ202;
    vdouble mv_detJ203;
    void setupWeight3d(const double w[], const uiint& numOfIntg, double w3d[]);
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, double Gzi[]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, double Gzi[]);
    void ShapeFunction8(vvdouble& N, const uiint& igauss,
                        const double& r, const double& s, const double& t);
    void ShapeFunction20(vvdouble& N, const uiint& igauss,
                         const double& r, const double& s, const double& t,
                         const double& r2, const double& s2, const double& t2);
    void ShapeDeriv8(vvvdouble& dNdr, const uiint& igauss,
                     const double& r, const double& s, const double& t);
    void ShapeDeriv20(vvvdouble& dNdr, const uiint& igauss,
                      const double& r, const double& s, const double& t,
                      const double& r2, const double& s2, const double& t2);
    vdouble mvIntegValue20;
    vdouble mvIntegValue8;
    void setupIntegValue20();
    void setupIntegValue8();
public:
    double& N81(const uiint& igauss, const uiint& ishape);
    double& N82(const uiint& igauss, const uiint& ishape);
    double& N201(const uiint& igauss, const uiint& ishape);
    double& N202(const uiint& igauss, const uiint& ishape);
    double& N203(const uiint& igauss, const uiint& ishape);
    vdouble& N81(const uiint& igauss);
    vdouble& N82(const uiint& igauss);
    vdouble& N201(const uiint& igauss);
    vdouble& N202(const uiint& igauss);
    vdouble& N203(const uiint& igauss);
    vvdouble& N81() {
        return mvN81;
    }
    vvdouble& N82() {
        return mvN82;
    }
    vvdouble& N201() {
        return mvN201;
    }
    vvdouble& N202() {
        return mvN202;
    }
    vvdouble& N203() {
        return mvN203;
    }
    double& dNdr81(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr82(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr201(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr202(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr203(const uiint& igauss, const uiint& ishape, const uiint& axis);
    vvdouble& dNdr81(const uiint& igauss);
    vvdouble& dNdr82(const uiint& igauss);
    vvdouble& dNdr201(const uiint& igauss);
    vvdouble& dNdr202(const uiint& igauss);
    vvdouble& dNdr203(const uiint& igauss);
    vvvdouble& dNdr81() {
        return mvdNdr81;
    }
    vvvdouble& dNdr82() {
        return mvdNdr82;
    }
    vvvdouble& dNdr201() {
        return mvdNdr201;
    }
    vvvdouble& dNdr202() {
        return mvdNdr202;
    }
    vvvdouble& dNdr203() {
        return mvdNdr203;
    }
    double* Weight3d(const uiint& integNum1D);
    double& Weight3dpt1();
    double& Weight3dpt2(const uiint& igauss);
    double& Weight3dpt3(const uiint& igauss);
    double& Weight3dpt4(const uiint& igauss);
    double& Weight3dpt5(const uiint& igauss);
    void Calc_dNdx8(const uiint& numOfInteg, CElement *pElement);
    void Calc_dNdx20(const uiint& numOfInteg, CElement *pElement);
    double& dNdx81(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return  mvdNdx81[igauss][ishape][axis];
    }
    double& dNdx82(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return  mvdNdx82[igauss][ishape][axis];
    }
    double& dNdx201(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx201[igauss][ishape][axis];
    }
    double& dNdx202(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx202[igauss][ishape][axis];
    }
    double& dNdx203(const uiint& igauss, const uiint& ishape, const uiint& axis) {
        return mvdNdx203[igauss][ishape][axis];
    }
    vvdouble& dNdx81(const uiint& igauss) {
        return  mvdNdx81[igauss];
    }
    vvdouble& dNdx82(const uiint& igauss) {
        return  mvdNdx82[igauss];
    }
    vvdouble& dNdx201(const uiint& igauss) {
        return mvdNdx201[igauss];
    }
    vvdouble& dNdx202(const uiint& igauss) {
        return mvdNdx202[igauss];
    }
    vvdouble& dNdx203(const uiint& igauss) {
        return mvdNdx203[igauss];
    }
    vvvdouble& dNdx81() {
        return  mvdNdx81;
    }
    vvvdouble& dNdx82() {
        return  mvdNdx82;
    }
    vvvdouble& dNdx201() {
        return mvdNdx201;
    }
    vvvdouble& dNdx202() {
        return mvdNdx202;
    }
    vvvdouble& dNdx203() {
        return mvdNdx203;
    }
    double& detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, const uiint& ishape);
    double& detJ81(const uiint& igauss) {
        return mv_detJ81[igauss];
    }
    double& detJ82(const uiint& igauss) {
        return mv_detJ82[igauss];
    }
    double& detJ201(const uiint& igauss) {
        return mv_detJ201[igauss];
    }
    double& detJ202(const uiint& igauss) {
        return mv_detJ202[igauss];
    }
    double& detJ203(const uiint& igauss) {
        return mv_detJ203[igauss];
    }
    vdouble& detJ81() {
        return  mv_detJ81;
    }
    vdouble& detJ82() {
        return  mv_detJ82;
    }
    vdouble& detJ201() {
        return mv_detJ201;
    }
    vdouble& detJ202() {
        return mv_detJ202;
    }
    vdouble& detJ203() {
        return mv_detJ203;
    }
    vdouble& getIntegralValue20() {
        return mvIntegValue20;
    }
    vdouble& getIntegralValue8() {
        return mvIntegValue8;
    }
    double& getIntegralValue20(const uiint& ishape) {
        return mvIntegValue20[ishape];
    }
    double& getIntegralValue8(const uiint& ishape) {
        return mvIntegValue8[ishape];
    }
};
#endif	/* _SHAPEHEXA_H */
}
