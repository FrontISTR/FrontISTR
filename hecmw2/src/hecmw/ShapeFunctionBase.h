/* 
 * File:   ShapeFunctionBase.h
 * Author: ktakeda
 *
 * Created on 2010/02/10, 15:26
 */
#include "TypeDef.h"

#include "Jacobian.h"

#include "ISTR2Edge.h"
#include "Edge2ISTR.h"

#include "Element.h"

namespace pmw{
#ifndef _SHAPEFUNCTIONBASE_H
#define	_SHAPEFUNCTIONBASE_H
class CShapeFunctionBase{
public:
    CShapeFunctionBase();
    virtual ~CShapeFunctionBase();
    
protected:
    //領域確保
    void ResizeShape(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape);
    void ResizeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof);
    void Resize2ndDeriv(v4double& d2Ndr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof);
    void Resize_detJ(vdouble& detJ, const uiint& numOfIntg);
    
    //積分点数の一覧
    vuint mvIntegNum;  //積分点の個数( Hexaの場合は,3次元に展開した積分点数 )
    vuint mvIntegNum1D;//Hexa,Quadの1次元方向の積分点の個数

    CISTR2Edge *mpISTR2Edge;
    CEdge2ISTR *mpEdge2ISTR;
    
    void ShapeDerivXYZ(const uiint& numOfInteg, const uiint& numOfShape, const uiint& dof,
                    const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx, vdouble& detJ);
public:
    // 積分点数
    uiint getNumOfIntegType(){ return mvIntegNum.size();}
    uiint& getNumOfInteg(const uiint& i){ return mvIntegNum[i];}
    // Hexaでの1次元方向の積分点数
    uiint getNumOfIntegType1D(){ return mvIntegNum1D.size();}
    uiint& getNumOfInteg1D(const uiint& i){ return mvIntegNum1D[i];}
};
#endif	/* _SHAPEFUNCTIONBASE_H */
}



