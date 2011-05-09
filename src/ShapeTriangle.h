/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeTriangle.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "ShapeFunctionBase.h"
namespace pmw{
#ifndef _SHAPETRIANGLE_H
#define	_SHAPETRIANGLE_H
class CShapeTriangle:public CShapeFunctionBase{
private:
    CShapeTriangle();
public:
    static CShapeTriangle* Instance(){
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
    void ShapeFunction3(vvdouble& N, const uint& igauss, const double& r, const double& s);
    void ShapeFunction6(vvdouble& N, const uint& igauss, const double& r, const double& s);
    void ShapeDeriv3(vvvdouble& dNdr, const uint& igauss);
    void ShapeDeriv6(vvvdouble& dNdr, const uint& igauss, const double& r, const double& s);
    void Shape_2ndDeriv3();
    void Shape_2ndDeriv6(v4double& d2Ndr, const uint& igauss);
public:
    double& N31(const uint& igauss, const uint& ishape);
    double& N63(const uint& igauss, const uint& ishape);
    vdouble& N31(const uint& igauss){ return mvN31[igauss];}
    vdouble& N63(const uint& igauss){ return mvN63[igauss];}
    vvdouble& N31(){ return mvN31;}
    vvdouble& N63(){ return mvN63;}
    double& dNdr31(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    double& dNdr63(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    vvdouble& dNdr31(const uint& igauss){ return mvdNdr31[igauss];}
    vvdouble& dNdr63(const uint& igauss){ return mvdNdr63[igauss];}
    vvvdouble& dNdr31(){ return mvdNdr31;}
    vvvdouble& dNdr63(){ return mvdNdr63;}
    double& d2Ndr31(const uint& igauss, const uint& ishape, const uint& iaxis, const uint& jaxis);
    double& d2Ndr63(const uint& igauss, const uint& ishape, const uint& iaxis, const uint& jaxis);
    vvvdouble& d2Ndr31(const uint& igauss){ return mvd2Ndr31[igauss];}
    vvvdouble& d2Ndr63(const uint& igauss){ return mvd2Ndr63[igauss];}
    v4double& d2Ndr31(){ return mvd2Ndr31;}
    v4double& d2Ndr63(){ return mvd2Ndr63;}
    double* Weight(const uint& integNum);
    double& Weight_pt1();
    double& Weight_pt3(const uint& igauss);
};
#endif	/* _SHAPETRIANGLE_H */
}
