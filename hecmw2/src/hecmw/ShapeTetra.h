/* 
 * File:   ShapeTetra.h
 * Author: ktakeda
 *
 * Created on 2010/01/28, 16:45
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
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[6];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[10];

    //積分座標(体積座標): L1,L2,L3
    // 
    double mGzi1[1][3];
    double mGzi4[4][3];
    double mGzi15[15][3];

    //積分点での重み
    double mW1[1];
    double mW4[4];
    double mW15[15];


    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN41;  //[1][4]    4節点要素  積分点1
    vvdouble mvN101; //[1][10]  10節点要素  積分点1
    vvdouble mvN104; //[4][10]  10節点要素  積分点4
    vvdouble mvN1015;//[15][10] 10節点要素 積分点15

    //導関数:dNdr[積分点][i][vol_coord]
    vvvdouble mvdNdr41;  //[1][4][3]    4節点要素  積分点1
    vvvdouble mvdNdr101; //[1][10][3]  10節点要素  積分点1
    vvvdouble mvdNdr104; //[4][10][3]  10節点要素  積分点4
    vvvdouble mvdNdr1015;//[15][10][3] 10節点要素 積分点15

    //ヤコビアン
    CJacobian *mpJacobi41;
    CJacobian *mpJacobi101;
    CJacobian *mpJacobi104;
    CJacobian *mpJacobi1015;
    //導関数(グローバル座標): dNdx[積分点][i][xyz]
    vvvdouble mvdNdx41; //[1][4][3]
    vvvdouble mvdNdx101;//[1][10][3]
    vvvdouble mvdNdx104;//[4][10][3]
    vvvdouble mvdNdx1015;//[15][10][3]

    //detJ : 配列[igauss]
    vdouble mv_detJ41;
    vdouble mv_detJ101;
    vdouble mv_detJ104;
    vdouble mv_detJ1015;


    
    //形状関数 & 導関数のセットアップ
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);



    //形状関数
    void ShapeFunction4(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& L3);
    void ShapeFunction10(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& L3);

    //導関数
    void ShapeDeriv4(vvvdouble& dNdr, const uiint& igauss);
    void ShapeDeriv10(vvvdouble& dNdr, const uiint& igauss,
                         const double& L1, const double& L2, const double& L3);

    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue10;
    vdouble mvIntegValue4;
    void setupIntegValue10();
    void setupIntegValue4();

public:
    // 形状関数   引数:積分点位置index, 形状関数番号(節点)
    double& N41(const uiint& igauss, const uiint& ishape);
    double& N101(const uiint& igauss, const uiint& ishape);
    double& N104(const uiint& igauss, const uiint& ishape);
    double& N1015(const uiint& igauss, const uiint& ishape);
    vdouble& N41(const uiint& igauss){ return mvN41[igauss];}
    vdouble& N101(const uiint& igauss){ return mvN101[igauss];}
    vdouble& N104(const uiint& igauss){ return mvN104[igauss];}
    vdouble& N1015(const uiint& igauss){ return mvN1015[igauss];}
    vvdouble& N41(){ return mvN41;}
    vvdouble& N101(){ return mvN101;}
    vvdouble& N104(){ return mvN104;}
    vvdouble& N1015(){ return mvN1015;}

    // 導関数   引数:積分点位置index, 形状関数番号(節点),微分方向(axis:0,1,2)
    double& dNdr41(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr101(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr104(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr1015(const uiint& igauss, const uiint& ishape, const uiint& axis);
    vvdouble& dNdr41(const uiint& igauss){ return mvdNdr41[igauss];}
    vvdouble& dNdr101(const uiint& igauss){ return mvdNdr101[igauss];}
    vvdouble& dNdr104(const uiint& igauss){ return mvdNdr104[igauss];}
    vvdouble& dNdr1015(const uiint& igauss){ return mvdNdr1015[igauss];}
    vvvdouble& dNdr41(){ return mvdNdr41;}
    vvvdouble& dNdr101(){ return mvdNdr101;}
    vvvdouble& dNdr104(){ return mvdNdr104;}
    vvvdouble& dNdr1015(){ return mvdNdr1015;}

    
    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uiint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt1();
    double& Weight_pt4(const uiint& igauss);
    double& Weight_pt15(const uiint& igauss);

    // dNdx 計算セットアップ
    void Calc_dNdx4(const uiint& numOfInteg, CElement *pElement);
    void Calc_dNdx10(const uiint& numOfInteg, CElement *pElement);

    // dNdx (空間座標の導関数)
    double& dNdx41(const uiint& igauss, const uiint& ishape, const uiint& axis){  return  mvdNdx41[igauss][ishape][axis];}
    double& dNdx101(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx101[igauss][ishape][axis];}
    double& dNdx104(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx104[igauss][ishape][axis];}
    double& dNdx1015(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx1015[igauss][ishape][axis];}

    vvdouble& dNdx41(const uiint& igauss){  return  mvdNdx41[igauss];}
    vvdouble& dNdx101(const uiint& igauss){ return mvdNdx101[igauss];}
    vvdouble& dNdx104(const uiint& igauss){ return mvdNdx104[igauss];}
    vvdouble& dNdx1015(const uiint& igauss){ return mvdNdx1015[igauss];}

    vvvdouble& dNdx41(){   return   mvdNdx41;}
    vvvdouble& dNdx101(){  return  mvdNdx101;}
    vvvdouble& dNdx104(){  return  mvdNdx104;}
    vvvdouble& dNdx1015(){ return mvdNdx1015;}

    // detJ : [積分点][ishape]
    double& detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss);

    double& detJ41(const uiint& igauss){ return mv_detJ41[igauss];}
    double& detJ101(const uiint& igauss){ return mv_detJ101[igauss];}
    double& detJ104(const uiint& igauss){ return mv_detJ104[igauss];}
    double& detJ1015(const uiint& igauss){ return mv_detJ1015[igauss];}
    vdouble& detJ41(){   return   mv_detJ41;}
    vdouble& detJ101(){  return  mv_detJ101;}
    vdouble& detJ104(){  return  mv_detJ104;}
    vdouble& detJ1015(){ return mv_detJ1015;}


    // 形状関数積分値
    vdouble& getIntegValue10(){ return mvIntegValue10;}
    vdouble& getIntegValue4(){ return mvIntegValue4;}
    double& getIntegValue10(const uiint& ishape){ return mvIntegValue10[ishape];}
    double& getIntegValue4(const uiint& ishape){ return mvIntegValue4[ishape];}
};
#endif	/* _SHAPETETRA_H */
}


