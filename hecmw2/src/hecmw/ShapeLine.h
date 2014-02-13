/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeLine.h
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
#ifndef _SHAPELINE_H
#define	_SHAPELINE_H
class CShapeLine:public CShapeFunctionBase
{
private:
    CShapeLine();
public:
    static CShapeLine* Instance() {
        static CShapeLine moShapeLine;
        return &moShapeLine;
    }
    virtual ~CShapeLine();
private:
    double mGzi1[1][1];
    double mGzi2[2][1];
    double mW1[1];
    double mW2[2];
    vvdouble mvN21;
    vvdouble mvN32;
    vvvdouble mvdNdr21;
    vvvdouble mvdNdr32;
    v4double mvd2Ndr21;
    v4double mvd2Ndr32;
    void setupShapeFunction();
    void setupShapeDeriv();
    void setupShape2ndDeriv();
    void ShapeFunction2(vvdouble& N, const uiint& igauss, const double& r);
    void ShapeFunction3(vvdouble& N, const uiint& igauss, const double& r);
    void ShapeDeriv2();
    void ShapeDeriv3(vvvdouble& dNdr, const uiint& igauss, const double& r);
    void Shape_2ndDeriv2();
    void Shape_2ndDeriv3();
    vdouble mvIntegValue3;
    vdouble mvIntegValue2;
    void setupIntegValue3();
    void setupIntegValue2();
public:
    double& N21(const uiint& igauss, const uiint& ishape);
    double& N32(const uiint& igauss, const uiint& ishape);
    vdouble& N21(const uiint& igauss) {
        return mvN21[igauss];
    }
    vdouble& N32(const uiint& igauss) {
        return mvN32[igauss];
    }
    vvdouble& N21() {
        return mvN21;
    }
    vvdouble& N32() {
        return mvN32;
    }
    double& dNdr21(const uiint& igauss, const uiint& ishape);
    double& dNdr32(const uiint& igauss, const uiint& ishape);
    vvdouble& dNdr21(const uiint& igauss) {
        return mvdNdr21[igauss];
    }
    vvdouble& dNdr32(const uiint& igauss) {
        return mvdNdr32[igauss];
    }
    vvvdouble& dNdr21() {
        return mvdNdr21;
    }
    vvvdouble& dNdr32() {
        return mvdNdr32;
    }
    double& d2Ndr21(const uiint& igauss, const uiint& ishape);
    double& d2Ndr32(const uiint& igauss, const uiint& ishape);
    vvvdouble& d2Ndr21(const uiint& igauss) {
        return mvd2Ndr21[igauss];
    }
    vvvdouble& d2Ndr32(const uiint& igauss) {
        return mvd2Ndr32[igauss];
    }
    v4double& d2Ndr21() {
        return mvd2Ndr21;
    }
    v4double& d2Ndr32() {
        return mvd2Ndr32;
    }
    double* Weight(const uiint& integNum);
    double& Weight_pt1();
    double& Weight_pt2(const uiint& igauss);
    vdouble& getIntegValue3() {
        return mvIntegValue3;
    }
    vdouble& getIntegValue2() {
        return mvIntegValue2;
    }
    double& getIntegValue3(const uiint& ishape) {
        return mvIntegValue3[ishape];
    }
    double& getIntegValue2(const uiint& ishape) {
        return mvIntegValue2[ishape];
    }
};
#endif	/* _SHAPELINE_H */
}
