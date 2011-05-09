/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ShapeFunctionBase.h
|
|                     Written by T.Takeda,    2010/06/01
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
namespace pmw{
#ifndef _SHAPEFUNCTIONBASE_H
#define	_SHAPEFUNCTIONBASE_H
class CShapeFunctionBase{
public:
    CShapeFunctionBase();
    virtual ~CShapeFunctionBase();
protected:
    void ResizeShape(vvdouble& N, const uint& numOfIntg, const uint& numOfShape);
    void ResizeDeriv(vvvdouble& dNdr, const uint& numOfIntg, const uint& numOfShape, const uint& dof);
    void Resize2ndDeriv(v4double& d2Ndr, const uint& numOfIntg, const uint& numOfShape, const uint& dof);
    void Resize_detJ(vdouble& detJ, const uint& numOfIntg);
    vuint mvIntegNum;  
    vuint mvIntegNum1D;
    CISTR2Edge *mpISTR2Edge;
    CEdge2ISTR *mpEdge2ISTR;
    void ShapeDerivXYZ(const uint& numOfInteg, const uint& numOfShape, const uint& dof,
                    const vvvdouble& dNdr, CJacobian* pJacobi, vector<CNode*>& vNode,vvvdouble& dNdx, vdouble& detJ);
public:
    uint getNumOfIntegType(){ return mvIntegNum.size();}
    uint& getNumOfInteg(const uint& i){ return mvIntegNum[i];}
    uint getNumOfIntegType1D(){ return mvIntegNum1D.size();}
    uint& getNumOfInteg1D(const uint& i){ return mvIntegNum1D[i];}
};
#endif	/* _SHAPEFUNCTIONBASE_H */
}
