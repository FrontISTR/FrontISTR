//
//  ShapeFunctionBase.cpp
//
//
//
//
//              2010.02.10
//              k.Takeda
#include "ShapeFunctionBase.h"
using namespace pmw;


CShapeFunctionBase::CShapeFunctionBase()
{
    mpISTR2Edge= CISTR2Edge::Instance();
    mpEdge2ISTR= CEdge2ISTR::Instance();
}
CShapeFunctionBase::~CShapeFunctionBase()
{
    ;
}

// ----
// 形状関数 領域確保
// ----
void CShapeFunctionBase::ResizeShape(vvdouble& N, const uint& numOfIntg, const uint& numOfShape)
{
    uint igauss;

    N.resize(numOfIntg);
    for(igauss=0; igauss< numOfIntg; igauss++){
        N[igauss].resize(numOfShape);
    };
}
// ----
// 導関数 領域確保
// ----
void CShapeFunctionBase::ResizeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const uint& dof)
{
    uint igauss,ishape;

    dNdr.resize(numOfIntg);
    for(igauss=0; igauss< numOfIntg; igauss++){
        dNdr[igauss].resize(numOfShape);
        for(ishape=0; ishape< numOfShape; ishape++){
            dNdr[igauss][ishape].resize(dof);
        };
    };
}
// ----
// 2次導関数 領域確保
// ----
void CShapeFunctionBase::Resize2ndDeriv(v4double& d2Ndr, const uint& numOfIntg, const uint& numOfShape, const uint& dof)
{
    uint igauss,ishape,ir;

    d2Ndr.resize(numOfIntg);
    for(igauss=0; igauss < numOfIntg; igauss++){
        d2Ndr[igauss].resize(numOfShape);

        for(ishape=0; ishape< numOfShape; ishape++){
            d2Ndr[igauss][ishape].resize(dof);

            for(ir=0; ir< dof; ir++){
                d2Ndr[igauss][ishape][ir].resize(dof);
            };
        };
    };
}
// ----
// detJ 配列の領域確保
// ----
void CShapeFunctionBase::Resize_detJ(vdouble& detJ, const uint& numOfIntg)
{
    uint igauss;

    detJ.resize(numOfIntg);
    
}

//  dNdxの計算 3次元(3*3):dof==3
//            2次元(2*2):dof==2
void CShapeFunctionBase::ShapeDerivXYZ(const uint& numOfInteg, const uint& numOfShape, const uint& dof,
                    const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx, vdouble& detJ)
{
    uint igauss,ishape,row,col;
    double val;//行列の加算結果 作業変数
    
    pJacobi->Calculate_J_invJ(dNdr, vNode);// J, invJ, detJ の計算

    for(igauss=0; igauss< numOfInteg; igauss++){
    for(ishape=0; ishape< numOfShape; ishape++){
        // dNdx,dNdy,dNdzの計算
        for(row=0; row< dof; row++){
            val=0.0;
            for(col=0; col< dof; col++){
                val += dNdr[igauss][ishape][col] * pJacobi->inverse_J(igauss,col,row);
            };
            dNdx[igauss][ishape][row]= val;
        };
        detJ[igauss]= pJacobi->detJ(igauss);
    };
    };
}









