/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeHexaNic.h
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
#ifndef _SHAPEHEXANIC_H
#define	_SHAPEHEXANIC_H
class CShapeHexaNic:public CShapeFunctionBase{
private:
    CShapeHexaNic();
public:
    static CShapeHexaNic* Instance(){
        static CShapeHexaNic moShapeHexaNic;
        return &moShapeHexaNic;
    }
    virtual ~CShapeHexaNic();
private:
    double mGzi1[1][3];
    double mGzi8[8][3];
    double mGzi27[27][3];
    double mW1[1];
    double mW8[8];
    double mW27[27];
    vvdouble mvN111; 
    vvdouble mvN118; 
    vvdouble mvN1127;
    vvvdouble mvdNdr111; 
    vvvdouble mvdNdr118; 
    vvvdouble mvdNdr1127;
    CJacobian *mpJacobi111;
    CJacobian *mpJacobi118;
    CJacobian *mpJacobi1127;
    vvvdouble mvdNdx111; 
    vvvdouble mvdNdx118; 
    vvvdouble mvdNdx1127;
    vdouble mv_detJ111;
    vdouble mv_detJ118;
    vdouble mv_detJ1127;
    void setupShapeFunction();
    void setupShapeDeriv();
    void ShapeFunction11(vvdouble& N, const uint& igauss, const double& r, const double& s, const double& t);
    void ShapeDeriv11(vvvdouble& dNdr, const uint& igauss, const double& r, const double& s, const double& t);
public:
    double& N111(const uint& igauss, const uint& ishape);
    double& N118(const uint& igauss, const uint& ishape);
    double& N1127(const uint& igauss, const uint& ishape);
    vdouble& N111(const uint& igauss);
    vdouble& N118(const uint& igauss);
    vdouble& N1127(const uint& igauss);
    vvdouble& N111(){ return mvN111;}
    vvdouble& N118(){ return mvN118;}
    vvdouble& N1127(){ return mvN1127;}
    double& dNdr111(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr118(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr1127(const uint& igauss, const uint& ishape, const uint& axis);
    vvdouble& dNdr111(const uint& igauss);
    vvdouble& dNdr118(const uint& igauss);
    vvdouble& dNdr1127(const uint& igauss);
    vvvdouble& dNdr111(){ return mvdNdr111;}
    vvvdouble& dNdr118(){ return mvdNdr118;}
    vvvdouble& dNdr1127(){ return mvdNdr1127;}
    double* Weight(const uint& integNum);
    double& Weight_pt1();
    double& Weight_pt8(const uint& igauss);
    double& Weight_pt27(const uint& igauss);
    void Calc_dNdx11(const uint& numOfInteg, CElement *pElement);
    double& dNdx111(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx111[igauss][ishape][axis];}
    double& dNdx118(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx118[igauss][ishape][axis];}
    double& dNdx1127(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx1127[igauss][ishape][axis];}
    vvdouble& dNdx111(const uint& igauss){  return  mvdNdx111[igauss];}
    vvdouble& dNdx118(const uint& igauss){  return  mvdNdx118[igauss];}
    vvdouble& dNdx1127(const uint& igauss){ return mvdNdx1127[igauss];}
    vvvdouble& dNdx111(){  return  mvdNdx111;}
    vvvdouble& dNdx118(){  return  mvdNdx118;}
    vvvdouble& dNdx1127(){ return mvdNdx1127;}
    double& detJ(const uint& elemType, const uint& numOfInteg, const uint& igauss);
    double& detJ111(const uint& igauss){ return mv_detJ111[igauss];}
    double& detJ118(const uint& igauss){ return mv_detJ118[igauss];}
    double& detJ1127(const uint& igauss){ return mv_detJ1127[igauss];}
    vdouble& detJ111(){  return  mv_detJ111;}
    vdouble& detJ118(){  return  mv_detJ118;}
    vdouble& detJ1127(){ return mv_detJ1127;}
};
#endif	/* _SHAPEHEXANIC_H */
}
