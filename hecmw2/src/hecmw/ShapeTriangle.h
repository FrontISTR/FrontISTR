/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeTriangle.h
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
#ifndef _SHAPETRIANGLE_H
#define	_SHAPETRIANGLE_H
class CShapeTriangle:public CShapeFunctionBase
{
private:
    CShapeTriangle();
public:
    static CShapeTriangle* Instance() {
        static CShapeTriangle moShapeTri;
        return &moShapeTri;
    }
    virtual ~CShapeTriangle();
private:
    double mGzi1[1][2];
    double mGzi3[3][2];
    double mW1[1];
    double mW3[3];
    vvdouble mvN31;
    vvdouble mvN63;
    vvvdouble mvdNdr31;
    vvvdouble mvdNdr63;
    v4double mvd2Ndr31;
    v4double mvd2Ndr63;
    void setupShapeFunction();
    void setupShapeDeriv();
    void setupShape2ndDeriv();
    void ShapeFunction3(vvdouble& N, const uiint& igauss, const double& r, const double& s);
    void ShapeFunction6(vvdouble& N, const uiint& igauss, const double& r, const double& s);
    void ShapeDeriv3(vvvdouble& dNdr, const uiint& igauss);
    void ShapeDeriv6(vvvdouble& dNdr, const uiint& igauss, const double& r, const double& s);
    void Shape_2ndDeriv3();
    void Shape_2ndDeriv6(v4double& d2Ndr, const uiint& igauss);
    vdouble mvIntegValue6;
    vdouble mvIntegValue3;
    void setupIntegValue6();
    void setupIntegValue3();
public:
    double& N31(const uiint& igauss, const uiint& ishape);
    double& N63(const uiint& igauss, const uiint& ishape);
    vdouble& N31(const uiint& igauss) {
        return mvN31[igauss];
    }
    vdouble& N63(const uiint& igauss) {
        return mvN63[igauss];
    }
    vvdouble& N31() {
        return mvN31;
    }
    vvdouble& N63() {
        return mvN63;
    }
    double& dNdr31(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    double& dNdr63(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    vvdouble& dNdr31(const uiint& igauss) {
        return mvdNdr31[igauss];
    }
    vvdouble& dNdr63(const uiint& igauss) {
        return mvdNdr63[igauss];
    }
    vvvdouble& dNdr31() {
        return mvdNdr31;
    }
    vvvdouble& dNdr63() {
        return mvdNdr63;
    }
    double& d2Ndr31(const uiint& igauss, const uiint& ishape, const uiint& iaxis, const uiint& jaxis);
    double& d2Ndr63(const uiint& igauss, const uiint& ishape, const uiint& iaxis, const uiint& jaxis);
    vvvdouble& d2Ndr31(const uiint& igauss) {
        return mvd2Ndr31[igauss];
    }
    vvvdouble& d2Ndr63(const uiint& igauss) {
        return mvd2Ndr63[igauss];
    }
    v4double& d2Ndr31() {
        return mvd2Ndr31;
    }
    v4double& d2Ndr63() {
        return mvd2Ndr63;
    }
    double* Weight(const uiint& integNum);
    double& Weight_pt1();
    double& Weight_pt3(const uiint& igauss);
    vdouble& getIntegValue6() {
        return mvIntegValue6;
    }
    vdouble& getIntegValue3() {
        return mvIntegValue3;
    }
    double& getIntegValue6(const uiint& ishape) {
        return mvIntegValue6[ishape];
    }
    double& getIntegValue3(const uiint& ishape) {
        return mvIntegValue3[ishape];
    }
};
#endif	/* _SHAPETRIANGLE_H */
}
