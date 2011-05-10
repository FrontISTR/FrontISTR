/* 
 * File:   ShapeQuad.h
 * Author: ktakeda
 *
 * Created on 2010/02/02, 16:15
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
//    // MW3 辺番号 -> FrontISTR 2次節点番号
//    uint mvEdge_2_ISTR[4];
//    // FrontISTR -> MW3 辺番号
//    uint mvISTR_2_Edge[8];

    //積分座標: gzi,eta
    double mGzi1[1][2];  //1点 積分
    double mGzi4[4][2];  //4点 積分
    double mGzi9[9][2];  //9点 積分

    //積分点の重み
    double mW1[1];  //1点 積分
    double mW4[4];  //4点 積分
    double mW9[9];  //9点 積分


    //形状関数:N[積分点][i] : i番の形状関数
    vvdouble mvN41;//[1][4]  4節点要素 積分点1
    vvdouble mvN84;//[4][8]  8節点要素 積分点4
    vvdouble mvN89;//[9][8]  8節点要素 積分点9

    //導関数: dNdr[積分点][i][coord]
    vvvdouble mvdNdr41;//[1][4][2]  4節点要素 積分点1
    vvvdouble mvdNdr84;//[4][8][2]  8節点要素 積分点4
    vvvdouble mvdNdr89;//[9][8][2]  8節点要素 積分点9
    //2次導関数: dNdr[積分点][i][coord][coord]
    v4double mvd2Ndr41;
    v4double mvd2Ndr84;
    v4double mvd2Ndr89;

    
    //形状関数 & 導関数のセットアップ
    void setupShapeFunction(vvdouble& N, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);
    void setupShapeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);
    void setupShape2ndDeriv(v4double& d2Ndr, const uint& numOfIntg, const uint& numOfShape, const double Gzi[][2]);


    //形状関数
    void ShapeFunction4(vvdouble& N, const uint& igauss,
                         const double& r, const double& s);
    void ShapeFunction8(vvdouble& N, const uint& igauss,
                         const double& r, const double& s);
    //導関数
    void ShapeDeriv4(vvvdouble& dNdr, const uint& igauss,
                        const double& r, const double& s);
    void ShapeDeriv8(vvvdouble& dNdr, const uint& igauss,
                         const double& r, const double& s);
    //2次導関数
    void Shape_2ndDeriv4();
    void Shape_2ndDeriv8(v4double& d2Ndr, const uint& igauss,
                                const double& r, const double& s);


    // Equivalent Node Force(等価節点力)
    //                  のための形状関数積分
    vdouble mvIntegValue8;
    vdouble mvIntegValue4;
    void setupIntegValue8();
    void setupIntegValue4();

public:
    // 形状関数   引数:積分点, 形状関数
    double& N41(const uint& igauss, const uint& ishape);
    double& N84(const uint& igauss, const uint& ishape);
    double& N89(const uint& igauss, const uint& ishape);
    vdouble& N41(const uint& igauss){ return mvN41[igauss];}
    vdouble& N84(const uint& igauss){ return mvN84[igauss];}
    vdouble& N89(const uint& igauss){ return mvN89[igauss];}
    vvdouble& N41(){ return mvN41;}
    vvdouble& N84(){ return mvN84;}
    vvdouble& N89(){ return mvN89;}
    // 導関数   引数:積分点, 形状関数,座標方向
    double& dNdr41(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    double& dNdr84(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    double& dNdr89(const uint& igauss, const uint& ishape, const uint& deriv_axis);
    vvdouble& dNdr41(const uint& igauss){ return mvdNdr41[igauss];}
    vvdouble& dNdr84(const uint& igauss){ return mvdNdr84[igauss];}
    vvdouble& dNdr89(const uint& igauss){ return mvdNdr89[igauss];}
    vvvdouble& dNdr41(){ return mvdNdr41;}
    vvvdouble& dNdr84(){ return mvdNdr84;}
    vvvdouble& dNdr89(){ return mvdNdr89;}
    // 2次導関数 引数:積分点, 形状関数, 座標方向, 座標方向
    double& d2Ndr41(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    double& d2Ndr84(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    double& d2Ndr89(const uint& igauss, const uint& ishape, const uint& deriv_axis0, const uint& deriv_axis1);
    vvvdouble& d2Ndr41(const uint& igauss){ return mvd2Ndr41[igauss];}
    vvvdouble& d2Ndr84(const uint& igauss){ return mvd2Ndr84[igauss];}
    vvvdouble& d2Ndr89(const uint& igauss){ return mvd2Ndr89[igauss];}
    v4double& d2Ndr41(){ return mvd2Ndr41;}
    v4double& d2Ndr84(){ return mvd2Ndr84;}
    v4double& d2Ndr89(){ return mvd2Ndr89;}


    // 返り値:重みの配列, 引数:積分点数
    double* Weight(const uint& integNum);
    // 返り値:重み, 引数:重みの配列Index
    double& Weight_pt1();
    double& Weight_pt4(const uint& igauss);
    double& Weight_pt9(const uint& igauss);


    // 形状関数の積分値
    vdouble& getIntegValue8(){ return mvIntegValue8;}
    vdouble& getIntegValue4(){ return mvIntegValue4;}
    double& getIntegValue8(const uint& ishape){ return mvIntegValue8[ishape];}
    double& getIntegValue4(const uint& ishape){ return mvIntegValue4[ishape];}
};
#endif	/* _SHAPEQUAD_H */
}











