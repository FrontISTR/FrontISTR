/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ShapeFunctionBase.h
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
#include "TypeDef.h"
#include "Jacobian.h"
#include "ISTR2Edge.h"
#include "Edge2ISTR.h"
#include "Element.h"
namespace pmw
{
#ifndef _SHAPEFUNCTIONBASE_H
#define	_SHAPEFUNCTIONBASE_H
class CShapeFunctionBase
{
public:
    CShapeFunctionBase();
    virtual ~CShapeFunctionBase();
protected:
    void ResizeShape(vvdouble& N, const uiint& numOfIntg, const uiint& numOfShape);
    void ResizeDeriv(vvvdouble& dNdr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof);
    void Resize2ndDeriv(v4double& d2Ndr, const uiint& numOfIntg, const uiint& numOfShape, const uiint& dof);
    void Resize_detJ(vdouble& detJ, const uiint& numOfIntg);
    vuint mvIntegNum;
    vuint mvIntegNum1D;
    CISTR2Edge *mpISTR2Edge;
    CEdge2ISTR *mpEdge2ISTR;
    void ShapeDerivXYZ(const uiint& numOfInteg, const uiint& numOfShape, const uiint& dof,
                       const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx, vdouble& detJ);
public:
    uiint getNumOfIntegType() {
        return mvIntegNum.size();
    }
    uiint& getNumOfInteg(const uiint& i) {
        return mvIntegNum[i];
    }
    uiint getNumOfIntegType1D() {
        return mvIntegNum1D.size();
    }
    uiint& getNumOfInteg1D(const uiint& i) {
        return mvIntegNum1D[i];
    }
};
#endif	/* _SHAPEFUNCTIONBASE_H */
}
