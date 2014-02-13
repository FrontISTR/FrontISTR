/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeFunctionBase.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
void CShapeFunctionBase::ResizeShape(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape)
{
    uiint igauss;
    N.resize(numOfIntg);
    for(igauss=0; igauss< numOfIntg; igauss++) {
        N[igauss].resize(numOfShape);
    };
}
void CShapeFunctionBase::ResizeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof)
{
    uiint igauss,ishape;
    dNdr.resize(numOfIntg);
    for(igauss=0; igauss< numOfIntg; igauss++) {
        dNdr[igauss].resize(numOfShape);
        for(ishape=0; ishape< numOfShape; ishape++) {
            dNdr[igauss][ishape].resize(dof);
        };
    };
}
void CShapeFunctionBase::Resize2ndDeriv(v4double& d2Ndr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof)
{
    uiint igauss,ishape,ir;
    d2Ndr.resize(numOfIntg);
    for(igauss=0; igauss < numOfIntg; igauss++) {
        d2Ndr[igauss].resize(numOfShape);
        for(ishape=0; ishape< numOfShape; ishape++) {
            d2Ndr[igauss][ishape].resize(dof);
            for(ir=0; ir< dof; ir++) {
                d2Ndr[igauss][ishape][ir].resize(dof);
            };
        };
    };
}
void CShapeFunctionBase::Resize_detJ(vdouble& detJ, const uiint& numOfIntg)
{
    uiint igauss;
    detJ.resize(numOfIntg);
}
void CShapeFunctionBase::ShapeDerivXYZ(const uiint& numOfInteg, const uiint& numOfShape, const uiint& dof,
                                       const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx, vdouble& detJ)
{
    uiint igauss,ishape,row,col;
    double val;
    pJacobi->Calculate_J_invJ(dNdr, vNode);
    for(igauss=0; igauss< numOfInteg; igauss++) {
        for(ishape=0; ishape< numOfShape; ishape++) {
            for(row=0; row< dof; row++) {
                val=0.0;
                for(col=0; col< dof; col++) {
                    val += dNdr[igauss][ishape][col] * pJacobi->inverse_J(igauss,col,row);
                };
                dNdx[igauss][ishape][row]= val;
            };
            detJ[igauss]= pJacobi->detJ(igauss);
        };
    };
}
