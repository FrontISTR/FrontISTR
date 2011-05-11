/* 
 * File:   ShapeLine.h
 * Author: ktakeda
 *
 * Created on 2010/02/03, 17:51
 */
#include "TypeDef.h"
#include "ShapeFunctionBase.h"

namespace pmw{
#ifndef _SHAPELINE_H
#define	_SHAPELINE_H
class CShapeLine:public CShapeFunctionBase{
private:
    CShapeLine();
public:
    static CShapeLine* Instance(){
        static CShapeLine moShapeLine;
        return &moShapeLine;
    }
    virtual ~CShapeLine();

private:
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[1];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[3];

    //積分座標: gzi
    double mGzi1[1][1];  //1点 積分
    double mGzi2[2][1];  //2点 積分

    //積分点の重み
    double mW1[1];  //1点 積分
    double mW2[2];  //2点 積分


    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN21;//[1][2] 2節点要素 積分点1
    vvdouble mvN32;//[2][3] 3節点要素 積分点2

    //導関数: dNdr[積分点][i]
    vvvdouble mvdNdr21;//[1][2][1]2節点要素 積分点1
    vvvdouble mvdNdr32;//[2][3][1]3節点要素 積分点2

    //2次導関数: dNdr[積分点][i]
    v4double mvd2Ndr21;//[1][2][1][1]
    v4double mvd2Ndr32;//[2][3][1][1]


    //形状関数 & 導関数のセットアップ
    void setupShapeFunction();
    void setupShapeDeriv();
    void setupShape2ndDeriv();


    //形状関数
    void ShapeFunction2(vvdouble& N, const uiint& igauss, const double& r);
    void ShapeFunction3(vvdouble& N, const uiint& igauss, const double& r);

    //導関数
    void ShapeDeriv2();
    //void ShapeDeriv3(double dNdr[][3][1], const uint& igauss, const double& r);
    void ShapeDeriv3(vvvdouble& dNdr, const uiint& igauss, const double& r);

    //2次導関数
    void Shape_2ndDeriv2();
    void Shape_2ndDeriv3();


    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue3;
    vdouble mvIntegValue2;
    void setupIntegValue3();
    void setupIntegValue2();

public:
    // 形状関数  引数:積分点, 形状関数
    double& N21(const uiint& igauss, const uiint& ishape);
    double& N32(const uiint& igauss, const uiint& ishape);
    vdouble& N21(const uiint& igauss){ return mvN21[igauss];}
    vdouble& N32(const uiint& igauss){ return mvN32[igauss];}
    vvdouble& N21(){ return mvN21;}
    vvdouble& N32(){ return mvN32;}
    // 導関数   引数:積分点, 形状関数
    double& dNdr21(const uiint& igauss, const uiint& ishape);
    double& dNdr32(const uiint& igauss, const uiint& ishape);
    vvdouble& dNdr21(const uiint& igauss){ return mvdNdr21[igauss];}
    vvdouble& dNdr32(const uiint& igauss){ return mvdNdr32[igauss];}
    vvvdouble& dNdr21(){ return mvdNdr21;}
    vvvdouble& dNdr32(){ return mvdNdr32;}
    // 2次導関数 引数:積分点, 形状関数
    double& d2Ndr21(const uiint& igauss, const uiint& ishape);
    double& d2Ndr32(const uiint& igauss, const uiint& ishape);
    vvvdouble& d2Ndr21(const uiint& igauss){ return mvd2Ndr21[igauss];}
    vvvdouble& d2Ndr32(const uiint& igauss){ return mvd2Ndr32[igauss];}
    v4double& d2Ndr21(){ return mvd2Ndr21;}
    v4double& d2Ndr32(){ return mvd2Ndr32;}

    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uiint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt1();
    double& Weight_pt2(const uiint& igauss);


    // 形状関数の積分値
    vdouble& getIntegValue3(){ return mvIntegValue3;}
    vdouble& getIntegValue2(){ return mvIntegValue2;}
    double& getIntegValue3(const uiint& ishape){ return mvIntegValue3[ishape];}
    double& getIntegValue2(const uiint& ishape){ return mvIntegValue2[ishape];}
};
#endif	/* _SHAPELINE_H */
}













