/* 
 * File:   ShapeHexaNic.h
 * Author: ktakeda
 *
 * Created on 2010/02/04, 14:28
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
    //積分座標:gzi,eta,zeta, 3D上の位置
    //
    double mGzi1[1][3];
    double mGzi8[8][3];
    double mGzi27[27][3];

    //積分点での重み
    double mW1[1];
    double mW8[8];
    double mW27[27];


    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN111; //[1][11]  11節点要素  積分点1
    vvdouble mvN118; //[8][11]  11節点要素  積分点8
    vvdouble mvN1127;//[27][11] 11節点要素  積分点27

    //導関数:dNdr[積分点][i][vol_coord]
    vvvdouble mvdNdr111; //[1][11][3]  11節点要素  積分点1
    vvvdouble mvdNdr118; //[8][11][3]  11節点要素  積分点8
    vvvdouble mvdNdr1127;//[27][11][3] 11節点要素  積分点27

    //ヤコビアン
    CJacobian *mpJacobi111;
    CJacobian *mpJacobi118;
    CJacobian *mpJacobi1127;
    //導関数(グローバル座標): dNdx[積分点][i][xyz]
    vvvdouble mvdNdx111; //[1][11][3]
    vvvdouble mvdNdx118; //[8][11][3]
    vvvdouble mvdNdx1127;//[27][11][3]
    //detJ
    vdouble mv_detJ111;
    vdouble mv_detJ118;
    vdouble mv_detJ1127;



    //形状関数 & 導関数のセットアップ
    void setupShapeFunction();
    void setupShapeDeriv();

    //形状関数
    //void ShapeFunction11(double N[][11], const uint& igauss, const double& r, const double& s, const double& t);
    void ShapeFunction11(vvdouble& N, const uint& igauss, const double& r, const double& s, const double& t);
    //導関数
    //void ShapeDeriv11(double dNdr[][11][3], const uint& igauss, const double& r, const double& s, const double& t);
    void ShapeDeriv11(vvvdouble& dNdr, const uint& igauss, const double& r, const double& s, const double& t);


public:
    // 形状関数   引数:積分点位置index, 形状関数番号(節点)
    double& N111(const uint& igauss, const uint& ishape);
    double& N118(const uint& igauss, const uint& ishape);
    double& N1127(const uint& igauss, const uint& ishape);
    vdouble& N111(const uint& igauss);
    vdouble& N118(const uint& igauss);
    vdouble& N1127(const uint& igauss);
    vvdouble& N111(){ return mvN111;}
    vvdouble& N118(){ return mvN118;}
    vvdouble& N1127(){ return mvN1127;}
    // 導関数   引数:積分点位置index, 形状関数番号(節点),微分方向(axis:0,1,2)
    double& dNdr111(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr118(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr1127(const uint& igauss, const uint& ishape, const uint& axis);
    vvdouble& dNdr111(const uint& igauss);
    vvdouble& dNdr118(const uint& igauss);
    vvdouble& dNdr1127(const uint& igauss);
    vvvdouble& dNdr111(){ return mvdNdr111;}
    vvvdouble& dNdr118(){ return mvdNdr118;}
    vvvdouble& dNdr1127(){ return mvdNdr1127;}

    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt1();
    double& Weight_pt8(const uint& igauss);
    double& Weight_pt27(const uint& igauss);

    // dNdx 計算セットアップ
    void Calc_dNdx11(const uint& numOfInteg, CElement *pElement);

    // dNdx (空間座標の導関数)
    double& dNdx111(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx111[igauss][ishape][axis];}
    double& dNdx118(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx118[igauss][ishape][axis];}
    double& dNdx1127(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx1127[igauss][ishape][axis];}

    vvdouble& dNdx111(const uint& igauss){  return  mvdNdx111[igauss];}
    vvdouble& dNdx118(const uint& igauss){  return  mvdNdx118[igauss];}
    vvdouble& dNdx1127(const uint& igauss){ return mvdNdx1127[igauss];}

    vvvdouble& dNdx111(){  return  mvdNdx111;}
    vvvdouble& dNdx118(){  return  mvdNdx118;}
    vvvdouble& dNdx1127(){ return mvdNdx1127;}

    // detJ : [積分点]
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












