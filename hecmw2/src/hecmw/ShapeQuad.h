/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeQuad.h
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
#ifndef _SHAPEQUAD_H
#define	_SHAPEQUAD_H
class CShapeQuad:public CShapeFunctionBase
{
private:
    CShapeQuad();
public:
    static CShapeQuad* Instance() {
        static CShapeQuad moShapeQuad;
        return &moShapeQuad;
    }
    virtual ~CShapeQuad();
private:
    double mGzi1[1][2];
    double mGzi4[4][2];
    double mGzi9[9][2];
    double mW1[1];
    double mW4[4];
    double mW9[9];
    vvdouble mvN41;
    vvdouble mvN84;
    vvdouble mvN89;
    vvvdouble mvdNdr41;
    vvvdouble mvdNdr84;
    vvvdouble mvdNdr89;
    v4double mvd2Ndr41;
    v4double mvd2Ndr84;
    v4double mvd2Ndr89;
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2]);
    void setupShape2ndDeriv(v4double& d2Ndr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][2]);
    void ShapeFunction4(vvdouble& N, const uiint& igauss,
                        const double& r, const double& s);
    void ShapeFunction8(vvdouble& N, const uiint& igauss,
                        const double& r, const double& s);
    void ShapeDeriv4(vvvdouble& dNdr, const uiint& igauss,
                     const double& r, const double& s);
    void ShapeDeriv8(vvvdouble& dNdr, const uiint& igauss,
                     const double& r, const double& s);
    void Shape_2ndDeriv4();
    void Shape_2ndDeriv8(v4double& d2Ndr, const uiint& igauss,
                         const double& r, const double& s);
    vdouble mvIntegValue8;
    vdouble mvIntegValue4;
    void setupIntegValue8();
    void setupIntegValue4();
public:
    double& N41(const uiint& igauss, const uiint& ishape);
    double& N84(const uiint& igauss, const uiint& ishape);
    double& N89(const uiint& igauss, const uiint& ishape);
    vdouble& N41(const uiint& igauss) {
        return mvN41[igauss];
    }
    vdouble& N84(const uiint& igauss) {
        return mvN84[igauss];
    }
    vdouble& N89(const uiint& igauss) {
        return mvN89[igauss];
    }
    vvdouble& N41() {
        return mvN41;
    }
    vvdouble& N84() {
        return mvN84;
    }
    vvdouble& N89() {
        return mvN89;
    }
    double& dNdr41(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    double& dNdr84(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    double& dNdr89(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    vvdouble& dNdr41(const uiint& igauss) {
        return mvdNdr41[igauss];
    }
    vvdouble& dNdr84(const uiint& igauss) {
        return mvdNdr84[igauss];
    }
    vvdouble& dNdr89(const uiint& igauss) {
        return mvdNdr89[igauss];
    }
    vvvdouble& dNdr41() {
        return mvdNdr41;
    }
    vvvdouble& dNdr84() {
        return mvdNdr84;
    }
    vvvdouble& dNdr89() {
        return mvdNdr89;
    }
    double& d2Ndr41(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1);
    double& d2Ndr84(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1);
    double& d2Ndr89(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis0, const uiint& deriv_axis1);
    vvvdouble& d2Ndr41(const uiint& igauss) {
        return mvd2Ndr41[igauss];
    }
    vvvdouble& d2Ndr84(const uiint& igauss) {
        return mvd2Ndr84[igauss];
    }
    vvvdouble& d2Ndr89(const uiint& igauss) {
        return mvd2Ndr89[igauss];
    }
    v4double& d2Ndr41() {
        return mvd2Ndr41;
    }
    v4double& d2Ndr84() {
        return mvd2Ndr84;
    }
    v4double& d2Ndr89() {
        return mvd2Ndr89;
    }
    double* Weight(const uiint& integNum);
    double& Weight_pt1();
    double& Weight_pt4(const uiint& igauss);
    double& Weight_pt9(const uiint& igauss);
    vdouble& getIntegValue8() {
        return mvIntegValue8;
    }
    vdouble& getIntegValue4() {
        return mvIntegValue4;
    }
    double& getIntegValue8(const uiint& ishape) {
        return mvIntegValue8[ishape];
    }
    double& getIntegValue4(const uiint& ishape) {
        return mvIntegValue4[ishape];
    }
};
#endif	/* _SHAPEQUAD_H */
}
