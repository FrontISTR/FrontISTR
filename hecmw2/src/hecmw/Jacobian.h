/* 
 * File:   Jacobian.h
 * Author: ktakeda
 *
 * Jacobian(3*3)の計算ルーチン
 *
 * Created on 2010/02/04, 16:32
 */
#include "TypeDef.h"
#include "Node.h"

#include "Logger.h"

namespace pmw{
#ifndef _JACOBIAN_H
#define	_JACOBIAN_H
class CJacobian{
public:
    CJacobian();
    virtual ~CJacobian();
    
private:
    double   m_detJ;// 3×3行列の行列式の値

    vvvdouble mvJ;   // [積分点]ごとの,3×3行列
    vvvdouble mvInvJ;// [積分点]ごとの,3×3逆行列
    vdouble mv_detJ; // [積分点]ごとの,3×3行列式値

    vvdouble mvJ33;  // 3×3行列
    vvdouble mvInv33;// 3×3行列の逆行列

public:
    void setupRegion(const uiint& numOfIntegPoint, const uiint& numOfShape);//領域確保(積分点数-形状関数のJacobian)
    void clearRegion(const uiint& numOfIntegPoint, const uiint& numOfShape);

    // ヤコビアン行列-逆行列を要素単位で丸ごと計算
    // ----
    // [積分点]別の,   J[3][3]をクラスメンバーmvJに保存
    // [積分点]別の,invJ[3][3]をクラスメンバーmvInvJに保存
    // [積分点]別の,      detJをクラスメンバーmv_detJに保存
    // ----
    void Calculate_J_invJ(const vvvdouble& dNdr, vector<CNode*>& vLocalNode);

    // ----
    // [積分点]別の"J","invJ","detJ"の提供
    // ----
    double& J(const uiint& igauss, const uiint& row, const uiint& col){
        return mvJ[igauss][row][col];
    }
    double& inverse_J(const uiint& igauss, const uiint& row, const uiint& col){
        return mvInvJ[igauss][row][col];
    }
    double& detJ(const uiint& igauss){
        return mv_detJ[igauss];
    }


    // ------
    // 3×3行列
    // ------
    double& detJ33(const vvdouble& Jmat);    
    vvdouble& inverse33(const vvdouble& Jmat);
    
};
#endif	/* _JACOBIAN_H */
}








