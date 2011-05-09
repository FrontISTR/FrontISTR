/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Jacobian.h
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
    double   m_detJ;
    vvvdouble mvJ;   
    vvvdouble mvInvJ;
    vdouble mv_detJ; 
    vvdouble mvJ33;  
    vvdouble mvInv33;
public:
    void setupRegion(const uint& numOfIntegPoint, const uint& numOfShape);
    void clearRegion(const uint& numOfIntegPoint, const uint& numOfShape);
    void Calculate_J_invJ(const vvvdouble& dNdr, vector<CNode*>& vLocalNode);
    double& J(const uint& igauss, const uint& row, const uint& col){
        return mvJ[igauss][row][col];
    }
    double& inverse_J(const uint& igauss, const uint& row, const uint& col){
        return mvInvJ[igauss][row][col];
    }
    double& detJ(const uint& igauss){
        return mv_detJ[igauss];
    }
    double& detJ33(const vvdouble& Jmat);    
    vvdouble& inverse33(const vvdouble& Jmat);
};
#endif	/* _JACOBIAN_H */
}
