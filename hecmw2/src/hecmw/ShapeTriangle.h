/* 
 * File:   ShapeTriangle.h
 * Author: ktakeda
 *
 * Created on 2010/02/02, 20:06
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
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[3];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[6];

    //積分座標: gzi,eta
    double mGzi1[1][2];  //1点 積分
    double mGzi3[3][2];  //3点 積分

    //積分点の重み
    double mW1[1];  //1点 積分
    double mW3[3];  //3点 積分


    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN31;//[1][3] 3節点要素 積分点1
    vvdouble mvN63;//[3][6] 6節点要素 積分点3

    //導関数: dNdr[積分点][i][coord]
    vvvdouble mvdNdr31;//[1][3][2] 3節点要素 積分点1
    vvvdouble mvdNdr63;//[3][6][2] 6節点要素 積分点3
    
    //2次導関数: dNdr[積分点][i][coord][coord]
    v4double mvd2Ndr31;//[1][3][2][2]
    v4double mvd2Ndr63;//[3][6][2][2]


    //形状関数 & 導関数のセットアップ
    void setupShapeFunction();
    void setupShapeDeriv();
    void setupShape2ndDeriv();

    
    //形状関数
    //void ShapeFunction3(double N[][3], const uint& igauss, const double& r, const double& s);
    //void ShapeFunction6(double N[][6], const uint& igauss, const double& r, const double& s);
    void ShapeFunction3(vvdouble& N, const uiint& igauss, const double& r, const double& s);
    void ShapeFunction6(vvdouble& N, const uiint& igauss, const double& r, const double& s);

    //導関数
    //void ShapeDeriv3(double dNdr[][3][2], const uint& igauss);
    //void ShapeDeriv6(double dNdr[][6][2], const uint& igauss, const double& r, const double& s);
    void ShapeDeriv3(vvvdouble& dNdr, const uiint& igauss);
    void ShapeDeriv6(vvvdouble& dNdr, const uiint& igauss, const double& r, const double& s);

    //2次導関数
    void Shape_2ndDeriv3();
    //void Shape_2ndDeriv6(double d2Ndr[][6][2][2], const uint& igauss);
    void Shape_2ndDeriv6(v4double& d2Ndr, const uiint& igauss);


    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue6;
    vdouble mvIntegValue3;
    void setupIntegValue6();
    void setupIntegValue3();

public:
    // 形状関数   引数:積分点, 形状関数
    double& N31(const uiint& igauss, const uiint& ishape);
    double& N63(const uiint& igauss, const uiint& ishape);
    vdouble& N31(const uiint& igauss){ return mvN31[igauss];}
    vdouble& N63(const uiint& igauss){ return mvN63[igauss];}
    vvdouble& N31(){ return mvN31;}
    vvdouble& N63(){ return mvN63;}
    // 導関数   引数:積分点, 形状関数,座標方向
    double& dNdr31(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    double& dNdr63(const uiint& igauss, const uiint& ishape, const uiint& deriv_axis);
    vvdouble& dNdr31(const uiint& igauss){ return mvdNdr31[igauss];}
    vvdouble& dNdr63(const uiint& igauss){ return mvdNdr63[igauss];}
    vvvdouble& dNdr31(){ return mvdNdr31;}
    vvvdouble& dNdr63(){ return mvdNdr63;}
    // 2次導関数 引数:積分点, 形状関数, 座標方向, 座標方向
    double& d2Ndr31(const uiint& igauss, const uiint& ishape, const uiint& iaxis, const uiint& jaxis);
    double& d2Ndr63(const uiint& igauss, const uiint& ishape, const uiint& iaxis, const uiint& jaxis);
    vvvdouble& d2Ndr31(const uiint& igauss){ return mvd2Ndr31[igauss];}
    vvvdouble& d2Ndr63(const uiint& igauss){ return mvd2Ndr63[igauss];}
    v4double& d2Ndr31(){ return mvd2Ndr31;}
    v4double& d2Ndr63(){ return mvd2Ndr63;}

    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uiint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt1();
    double& Weight_pt3(const uiint& igauss);


    // 形状関数の積分値
    vdouble& getIntegValue6(){ return mvIntegValue6;}
    vdouble& getIntegValue3(){ return mvIntegValue3;}
    double& getIntegValue6(const uiint& ishape){ return mvIntegValue6[ishape];}
    double& getIntegValue3(const uiint& ishape){ return mvIntegValue3[ishape];}

};
#endif	/* _SHAPETRIANGLE_H */
}

