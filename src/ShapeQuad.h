/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeQuad.h
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
#ifndef _SHAPEQUAD_H
#define	_SHAPEQUAD_H
class CShapeQuad:public CShapeFunctionBase{
private:
    CShapeQuad();
public:
    static CShapeQuad* Instance(){
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
    void setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);
    void setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);
    void setupShape2ndDeriv(v4double& d2Ndr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);
    void ShapeFunction4(vvdouble& N, const uint& igauss,
                         const double& r, const double& s);
    void ShapeFunction8(vvdouble& N, const uint& igauss,
                         const double& r, const double& s);
    void ShapeDeriv4(vvvdouble& dNdr, const uint& igauss,
                        const double& r, const double& s);
    void ShapeDeriv8(vvvdouble& dNdr, const uint& igauss,
                         const double& r, const double& s);
    void Shape_2ndDeriv4();
    void Shape_2ndDeriv8(v4double& d2Ndr, const uint& igauss,
                                const double& r, const double& s);
public:
    double& N41(const uint& igauss, const uint& ishape);
    double& N84(const uint& igauss, const uint& ishape);
    double& N89(const uint& igauss, const uint& ishape);
    vdouble& N41(const uint& igauss){ return mvN41[igauss];}
    vdouble& N84(const uint& igauss){ return mvN84[igauss];}
    vdouble& N89(const uint& igauss){ return mvN89[igauss];}
    vvdouble& N41(){ return mvN41;}
    vvdouble& N84(){ return mvN84;}
    vvdouble& N89(){ return mvN89;}
    double& dNdr41(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    double& dNdr84(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    double& dNdr89(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    vvdouble& dNdr41(const uint& igauss){ return mvdNdr41[igauss];}
    vvdouble& dNdr84(const uint& igauss){ return mvdNdr84[igauss];}
    vvdouble& dNdr89(const uint& igauss){ return mvdNdr89[igauss];}
    vvvdouble& dNdr41(){ return mvdNdr41;}
    vvvdouble& dNdr84(){ return mvdNdr84;}
    vvvdouble& dNdr89(){ return mvdNdr89;}
    double& d2Ndr41(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    double& d2Ndr84(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    double& d2Ndr89(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    vvvdouble& d2Ndr41(const uint& igauss){ return mvd2Ndr41[igauss];}
    vvvdouble& d2Ndr84(const uint& igauss){ return mvd2Ndr84[igauss];}
    vvvdouble& d2Ndr89(const uint& igauss){ return mvd2Ndr89[igauss];}
    v4double& d2Ndr41(){ return mvd2Ndr41;}
    v4double& d2Ndr84(){ return mvd2Ndr84;}
    v4double& d2Ndr89(){ return mvd2Ndr89;}
    double* Weight(const uint& integNum);
    double& Weight_pt1();
    double& Weight_pt4(const uint& igauss);
    double& Weight_pt9(const uint& igauss);
};
#endif	/* _SHAPEQUAD_H */
}
