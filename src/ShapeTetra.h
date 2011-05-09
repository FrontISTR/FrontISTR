/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeTetra.h
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
#ifndef _SHAPETETRA_H
#define	_SHAPETETRA_H
class CShapeTetra:public CShapeFunctionBase{
private:
    CShapeTetra();
public:
    static CShapeTetra* Instance(){
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
    void setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][3]);
    void setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][3]);
    void ShapeFunction4(vvdouble& N, const uint& igauss,
                         const double& L1, const double& L2, const double& L3);
    void ShapeFunction10(vvdouble& N, const uint& igauss,
                         const double& L1, const double& L2, const double& L3);
    void ShapeDeriv4(vvvdouble& dNdr, const uint& igauss);
    void ShapeDeriv10(vvvdouble& dNdr, const uint& igauss,
                         const double& L1, const double& L2, const double& L3);
public:
    double& N41(const uint& igauss, const uint& ishape);
    double& N101(const uint& igauss, const uint& ishape);
    double& N104(const uint& igauss, const uint& ishape);
    double& N1015(const uint& igauss, const uint& ishape);
    vdouble& N41(const uint& igauss){ return mvN41[igauss];}
    vdouble& N101(const uint& igauss){ return mvN101[igauss];}
    vdouble& N104(const uint& igauss){ return mvN104[igauss];}
    vdouble& N1015(const uint& igauss){ return mvN1015[igauss];}
    vvdouble& N41(){ return mvN41;}
    vvdouble& N101(){ return mvN101;}
    vvdouble& N104(){ return mvN104;}
    vvdouble& N1015(){ return mvN1015;}
    double& dNdr41(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr101(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr104(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr1015(const uint& igauss, const uint& ishape, const uint& axis);
    vvdouble& dNdr41(const uint& igauss){ return mvdNdr41[igauss];}
    vvdouble& dNdr101(const uint& igauss){ return mvdNdr101[igauss];}
    vvdouble& dNdr104(const uint& igauss){ return mvdNdr104[igauss];}
    vvdouble& dNdr1015(const uint& igauss){ return mvdNdr1015[igauss];}
    vvvdouble& dNdr41(){ return mvdNdr41;}
    vvvdouble& dNdr101(){ return mvdNdr101;}
    vvvdouble& dNdr104(){ return mvdNdr104;}
    vvvdouble& dNdr1015(){ return mvdNdr1015;}
    double* Weight(const uint& integNum);
    double& Weight_pt1();
    double& Weight_pt4(const uint& igauss);
    double& Weight_pt15(const uint& igauss);
    void Calc_dNdx4(const uint& numOfInteg, CElement *pElement);
    void Calc_dNdx10(const uint& numOfInteg, CElement *pElement);
    double& dNdx41(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx41[igauss][ishape][axis];}
    double& dNdx101(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx101[igauss][ishape][axis];}
    double& dNdx104(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx104[igauss][ishape][axis];}
    double& dNdx1015(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx1015[igauss][ishape][axis];}
    vvdouble& dNdx41(const uint& igauss){  return  mvdNdx41[igauss];}
    vvdouble& dNdx101(const uint& igauss){ return mvdNdx101[igauss];}
    vvdouble& dNdx104(const uint& igauss){ return mvdNdx104[igauss];}
    vvdouble& dNdx1015(const uint& igauss){ return mvdNdx1015[igauss];}
    vvvdouble& dNdx41(){   return   mvdNdx41;}
    vvvdouble& dNdx101(){  return  mvdNdx101;}
    vvvdouble& dNdx104(){  return  mvdNdx104;}
    vvvdouble& dNdx1015(){ return mvdNdx1015;}
    double& detJ(const uint& elemType, const uint& numOfInteg, const uint& igauss);
    double& detJ41(const uint& igauss){ return mv_detJ41[igauss];}
    double& detJ101(const uint& igauss){ return mv_detJ101[igauss];}
    double& detJ104(const uint& igauss){ return mv_detJ104[igauss];}
    double& detJ1015(const uint& igauss){ return mv_detJ1015[igauss];}
    vdouble& detJ41(){   return   mv_detJ41;}
    vdouble& detJ101(){  return  mv_detJ101;}
    vdouble& detJ104(){  return  mv_detJ104;}
    vdouble& detJ1015(){ return mv_detJ1015;}
};
#endif	/* _SHAPETETRA_H */
}
