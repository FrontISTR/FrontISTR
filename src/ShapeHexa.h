/* 
 * File:   ShapeHexa.h
 * Author: ktakeda
 *
 * Hexaの形状関数クラス
 *
 * Created on 2010/01/21, 18:41
 */
#include "TypeDef.h"
#include "ShapeFunctionBase.h"


namespace pmw{
#ifndef _SHAPEHEXA_H
#define	_SHAPEHEXA_H

class CShapeHexa:public CShapeFunctionBase{
private:
    CShapeHexa();
public:
    virtual ~CShapeHexa();

    static CShapeHexa* Instance(){
        static CShapeHexa moShapeHexa;
        return &moShapeHexa;
    }

private:
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[12];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[20];

    //積分座標(自然座標,正規化座標):1次元で表記
    double mGzi1[1];
    double mGzi2[2];
    double mGzi3[3];
    double mGzi4[4];
    double mGzi5[5];

    //積分点での重み 1次元
    double mW1[1];
    double mW2[2];
    double mW3[3];
    double mW4[4];
    double mW5[5];
    //積分点での重み 3次元: w[ir]*w[is]*w[it]
    double mW3d1[1];
    double mW3d2[8];
    double mW3d3[27];
    double mW3d4[64];
    double mW3d5[125];


    //形状関数:N[積分点][i] : i番の形状関数
    //
    vvdouble mvN81; //[1][8]     8節点要素 積分点1
    vvdouble mvN82; //[8][8]     8節点要素 積分点2(8)
    vvdouble mvN201;//[1]20]    20節点要素 積分点1
    vvdouble mvN202;//[8][20]   20節点要素 積分点2(8)
    vvdouble mvN203;//[27][20]  20節点要素 積分点3(27)
    
    //導関数: dNdr[積分点][i][rst]
    vvvdouble mvdNdr81; //[1][8][3]     8節点要素 積分点1
    vvvdouble mvdNdr82; //[8][8][3]     8節点要素 積分点2(8)
    vvvdouble mvdNdr201;//[1][20][3]   20節点要素 積分点1
    vvvdouble mvdNdr202;//[8][20][3]   20節点要素 積分点2(8)
    vvvdouble mvdNdr203;//[27][20][3]  20節点要素 積分点3(27)

    //ヤコビアン
    CJacobian *mpJacobi81;
    CJacobian *mpJacobi82;
    CJacobian *mpJacobi201;
    CJacobian *mpJacobi202;
    CJacobian *mpJacobi203;
    //導関数(グローバル座標): dNdx[積分点][i][xyz]
    vvvdouble mvdNdx81; //[1][8][3]
    vvvdouble mvdNdx82; //[8][8][3]
    vvvdouble mvdNdx201;//[1][20][3]
    vvvdouble mvdNdx202;//[8][20][3]
    vvvdouble mvdNdx203;//[27][20][3]
    // detJ : [積分点][ishape]のdetJ
    vdouble mv_detJ81;
    vdouble mv_detJ82;
    vdouble mv_detJ201;
    vdouble mv_detJ202;
    vdouble mv_detJ203;


    //重みの3D化
    void setupWeight3d(const double w[], const uint& numOfIntg, double w3d[]);
    


    //形状関数 & 導関数のセットアップ(20節点のShapeFunction,ShapeDeriv呼び出し)
    void setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, double Gzi[]);
    void setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, double Gzi[]);
    
    //形状関数計算ルーチン
    void ShapeFunction8(vvdouble& N, const uint& igauss,
                         const double& r, const double& s, const double& t);
    void ShapeFunction20(vvdouble& N, const uint& igauss,
                         const double& r, const double& s, const double& t,
                         const double& r2, const double& s2, const double& t2);
    
    //dNdr計算ルーチン
    void ShapeDeriv8(vvvdouble& dNdr, const uint& igauss,
                         const double& r, const double& s, const double& t);
    void ShapeDeriv20(vvvdouble& dNdr, const uint& igauss,
                         const double& r, const double& s, const double& t,
                         const double& r2, const double& s2, const double& t2);

    // CShapeFunctionBaseに移動
    //
    ////dNdx計算ルーチン (Calc_dNdx**からコール)
    //void ShapeDerivXYZ(const uint& numOfInteg, const uint& numOfShape, const uint& dof,
    //        const vvvdouble& dNdr, CJacobian *pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx);



    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue20;
    vdouble mvIntegValue8;
    void setupIntegValue20();
    void setupIntegValue8();
    
public:
    // 形状関数   引数:積分点位置3次元上のindex, 形状関数番号(節点)
    double& N81(const uint& igauss, const uint& ishape);
    double& N82(const uint& igauss, const uint& ishape);
    double& N201(const uint& igauss, const uint& ishape);
    double& N202(const uint& igauss, const uint& ishape);
    double& N203(const uint& igauss, const uint& ishape);

    vdouble& N81(const uint& igauss);
    vdouble& N82(const uint& igauss);
    vdouble& N201(const uint& igauss);
    vdouble& N202(const uint& igauss);
    vdouble& N203(const uint& igauss);

    vvdouble& N81(){ return mvN81;}
    vvdouble& N82(){ return mvN82;}
    vvdouble& N201(){ return mvN201;}
    vvdouble& N202(){ return mvN202;}
    vvdouble& N203(){ return mvN203;}

    // 導関数   引数:積分点位置3次元上のindex, 形状関数番号(節点),微分方向(axis:0,1,2)
    double& dNdr81(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr82(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr201(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr202(const uint& igauss, const uint& ishape, const uint& axis);
    double& dNdr203(const uint& igauss, const uint& ishape, const uint& axis);

    vvdouble& dNdr81(const uint& igauss);
    vvdouble& dNdr82(const uint& igauss);
    vvdouble& dNdr201(const uint& igauss);
    vvdouble& dNdr202(const uint& igauss);
    vvdouble& dNdr203(const uint& igauss);

    vvvdouble& dNdr81(){ return mvdNdr81;}
    vvvdouble& dNdr82(){ return mvdNdr82;}
    vvvdouble& dNdr201(){ return mvdNdr201;}
    vvvdouble& dNdr202(){ return mvdNdr202;}
    vvvdouble& dNdr203(){ return mvdNdr203;}

    // 返り値:重みの配列, 引数:一次元での積分点数
    double* Weight3d(const uint& integNum1D);
    // 返り値:3D上のIndexでの重み, 引数:重みの配列Index
    double& Weight3dpt1();
    double& Weight3dpt2(const uint& igauss);
    double& Weight3dpt3(const uint& igauss);
    double& Weight3dpt4(const uint& igauss);
    double& Weight3dpt5(const uint& igauss);

    
    // dNdx 計算セットアップ
    void Calc_dNdx8(const uint& numOfInteg, CElement *pElement);
    void Calc_dNdx20(const uint& numOfInteg, CElement *pElement);
    //void Calc_dNdx8(const uint& numOfInteg,vector<CNode*>& vNode);
    //void Calc_dNdx20(const uint& numOfInteg, vector<CNode*>& vVertNode, vector<CNode*>& vEdgeNode);

    // dNdx (空間座標の導関数)
    double& dNdx81(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx81[igauss][ishape][axis];}
    double& dNdx82(const uint& igauss, const uint& ishape, const uint& axis){  return  mvdNdx82[igauss][ishape][axis];}
    double& dNdx201(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx201[igauss][ishape][axis];}
    double& dNdx202(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx202[igauss][ishape][axis];}
    double& dNdx203(const uint& igauss, const uint& ishape, const uint& axis){ return mvdNdx203[igauss][ishape][axis];}

    vvdouble& dNdx81(const uint& igauss){  return  mvdNdx81[igauss];}
    vvdouble& dNdx82(const uint& igauss){  return  mvdNdx82[igauss];}
    vvdouble& dNdx201(const uint& igauss){ return mvdNdx201[igauss];}
    vvdouble& dNdx202(const uint& igauss){ return mvdNdx202[igauss];}
    vvdouble& dNdx203(const uint& igauss){ return mvdNdx203[igauss];}

    vvvdouble& dNdx81(){  return  mvdNdx81;}
    vvvdouble& dNdx82(){  return  mvdNdx82;}
    vvvdouble& dNdx201(){ return mvdNdx201;}
    vvvdouble& dNdx202(){ return mvdNdx202;}
    vvvdouble& dNdx203(){ return mvdNdx203;}
    
    // detJ : [積分点][ishape]
    double& detJ(const uint& elemType, const uint& numOfInteg, const uint& igauss, const uint& ishape);

    double& detJ81(const uint& igauss){ return mv_detJ81[igauss];}
    double& detJ82(const uint& igauss){ return mv_detJ82[igauss];}
    double& detJ201(const uint& igauss){ return mv_detJ201[igauss];}
    double& detJ202(const uint& igauss){ return mv_detJ202[igauss];}
    double& detJ203(const uint& igauss){ return mv_detJ203[igauss];}
    vdouble& detJ81(){  return  mv_detJ81;}
    vdouble& detJ82(){  return  mv_detJ82;}
    vdouble& detJ201(){ return mv_detJ201;}
    vdouble& detJ202(){ return mv_detJ202;}
    vdouble& detJ203(){ return mv_detJ203;}


    // 形状関数積分値
    vdouble& getIntegralValue20(){ return mvIntegValue20;}
    vdouble& getIntegralValue8(){ return mvIntegValue8;}
    double& getIntegralValue20(const uint& ishape){ return mvIntegValue20[ishape];}
    double& getIntegralValue8(const uint& ishape){ return mvIntegValue8[ishape];}
};
#endif	/* _SHAPEHEXA_H */
}












