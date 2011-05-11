/* 
 * File:   ShapePrism.h
 * Author: ktakeda
 *
 * Created on 2010/01/29, 20:17
 */
#include "TypeDef.h"
#include "ShapeFunctionBase.h"

namespace pmw{
#ifndef _SHAPEPRISM_H
#define	_SHAPEPRISM_H
class CShapePrism:public CShapeFunctionBase{
private:
    CShapePrism();
public:
    static CShapePrism* Instance(){
        static CShapePrism moShapePrism;
        return &moShapePrism;
    }
    virtual ~CShapePrism();

private:
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[9];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[15];

    //積分座標(面積座標,gzi): L1,L2,gzi
    //
    double mGzi2[2][3];  //2点 積分
    double mGzi6[6][3];  //6点 積分
    double mGzi9[9][3];  //9点 積分
    double mGzi18[18][3];//18点 積分

    //積分点の重み
    double mW2[2];  //2点 積分
    double mW6[6];  //6点 積分
    double mW9[9];  //9点 積分
    double mW18[18];//18点 積分

    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN62;   //[2][6]    6節点要素 積分点2
    vvdouble mvN156;  //[6][15]  15節点要素 積分点6
    vvdouble mvN159;  //[9][15]  15節点要素 積分点9
    vvdouble mvN1518; //[18][15] 15節点要素 積分点18

    //導関数: dNdr[積分点][i][integ_coord]
    vvvdouble mvdNdr62;  //[2][6][3]    6節点要素 積分点2
    vvvdouble mvdNdr156; //[6][15][3]  15節点要素 積分点6
    vvvdouble mvdNdr159; //[9][15][3]  15節点要素 積分点9
    vvvdouble mvdNdr1518;//[18][15][3] 15節点要素 積分点18

    //ヤコビアン
    CJacobian *mpJacobi62;
    CJacobian *mpJacobi156;
    CJacobian *mpJacobi159;
    CJacobian *mpJacobi1518;
    //導関数(グローバル座標): dNdx[積分点][i][xyz]
    vvvdouble mvdNdx62; //[2][6][3]
    vvvdouble mvdNdx156;//[6][15][3]
    vvvdouble mvdNdx159;//[9][15][3]
    vvvdouble mvdNdx1518;//[18][15][3]

    // detJ : 配列[igauss][ishape]
    vdouble mv_detJ62;
    vdouble mv_detJ156;
    vdouble mv_detJ159;
    vdouble mv_detJ1518;


    
    //形状関数 & 導関数のセットアップ
    void setupShapeFunction(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);
    void setupShapeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const double Gzi[][3]);


    //形状関数
    void ShapeFunction6(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& gzi);
    void ShapeFunction15(vvdouble& N, const uiint& igauss,
                         const double& L1, const double& L2, const double& gzi);

    //導関数
    void ShapeDeriv6(vvvdouble& dNdr, const uiint& igauss,
                        const double& L1, const double& L2, const double& gzi);
    void ShapeDeriv15(vvvdouble& dNdr, const uiint& igauss,
                         const double& L1, const double& L2, const double& gzi);


    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue15;
    vdouble mvIntegValue6;
    void setupIntegValue15();
    void setupIntegValue6();


public:
    // 形状関数   引数:積分点位置index, 形状関数番号(節点)
    double& N62(const uiint& igauss, const uiint& ishape);
    double& N156(const uiint& igauss, const uiint& ishape);
    double& N159(const uiint& igauss, const uiint& ishape);
    double& N1518(const uiint& igauss, const uiint& ishape);
    vdouble& N62(const uiint& igauss){ return mvN62[igauss];}
    vdouble& N156(const uiint& igauss){ return mvN156[igauss];}
    vdouble& N159(const uiint& igauss){ return mvN159[igauss];}
    vdouble& N1518(const uiint& igauss){ return mvN1518[igauss];}
    vvdouble& N62(){ return mvN62;}
    vvdouble& N156(){ return mvN156;}
    vvdouble& N159(){ return mvN159;}
    vvdouble& N1518(){ return mvN1518;}

    // 導関数   引数:積分点位置index, 形状関数番号(節点),微分方向(axis:0,1,2)
    double& dNdr62(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr156(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr159(const uiint& igauss, const uiint& ishape, const uiint& axis);
    double& dNdr1518(const uiint& igauss, const uiint& ishape, const uiint& axis);
    vvdouble& dNdr62(const uiint& igauss){ return mvdNdr62[igauss];}
    vvdouble& dNdr156(const uiint& igauss){ return mvdNdr156[igauss];}
    vvdouble& dNdr159(const uiint& igauss){ return mvdNdr159[igauss];}
    vvdouble& dNdr1518(const uiint& igauss){ return mvdNdr1518[igauss];}
    vvvdouble& dNdr62(){ return mvdNdr62;}
    vvvdouble& dNdr156(){ return mvdNdr156;}
    vvvdouble& dNdr159(){ return mvdNdr159;}
    vvvdouble& dNdr1518(){ return mvdNdr1518;}

    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uiint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt2(const uiint& igauss);
    double& Weight_pt6(const uiint& igauss);
    double& Weight_pt9(const uiint& igauss);
    double& Weight_pt18(const uiint& igauss);

    // dNdx 計算セットアップ
    void Calc_dNdx6(const uiint& numOfInteg, CElement *pElement);
    void Calc_dNdx15(const uiint& numOfInteg, CElement *pElement);

    // dNdx (空間座標の導関数)
    double& dNdx62(const uiint& igauss, const uiint& ishape, const uiint& axis){  return  mvdNdx62[igauss][ishape][axis];}
    double& dNdx156(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx156[igauss][ishape][axis];}
    double& dNdx159(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx159[igauss][ishape][axis];}
    double& dNdx1518(const uiint& igauss, const uiint& ishape, const uiint& axis){ return mvdNdx1518[igauss][ishape][axis];}

    vvdouble& dNdx62(const uiint& igauss){  return  mvdNdx62[igauss];}
    vvdouble& dNdx156(const uiint& igauss){ return mvdNdx156[igauss];}
    vvdouble& dNdx159(const uiint& igauss){ return mvdNdx159[igauss];}
    vvdouble& dNdx1518(const uiint& igauss){ return mvdNdx1518[igauss];}

    vvvdouble& dNdx62(){   return   mvdNdx62;}
    vvvdouble& dNdx156(){  return  mvdNdx156;}
    vvvdouble& dNdx159(){  return  mvdNdx159;}
    vvvdouble& dNdx1518(){ return mvdNdx1518;}

    // detJ : [積分点][ishape]
    double& detJ(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss);

    double& detJ62(const uiint& igauss){ return mv_detJ62[igauss];}
    double& detJ156(const uiint& igauss){ return mv_detJ156[igauss];}
    double& detJ159(const uiint& igauss){ return mv_detJ159[igauss];}
    double& detJ1518(const uiint& igauss){ return mv_detJ1518[igauss];}
    vdouble& detJ62(){   return   mv_detJ62;}
    vdouble& detJ156(){  return  mv_detJ156;}
    vdouble& detJ159(){  return  mv_detJ159;}
    vdouble& detJ1518(){ return mv_detJ1518;}


    // 形状関数積分値
    vdouble& getIntegValue15(){ return mvIntegValue15;}
    vdouble& getIntegValue6(){ return mvIntegValue6;}
    double& getIntegValue15(const uiint& ishape){ return mvIntegValue15[ishape];}
    double& getIntegValue6(const uiint& ishape){ return mvIntegValue6[ishape];}
};
#endif	/* _SHAPEPRISM_H */
}



