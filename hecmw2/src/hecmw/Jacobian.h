/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Jacobian.h
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
#include "Node.h"
#include "Logger.h"
namespace pmw
{
#ifndef _JACOBIAN_H
#define	_JACOBIAN_H
class CJacobian
{
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
    void setupRegion(const uiint& numOfIntegPoint, const uiint& numOfShape);
    void clearRegion(const uiint& numOfIntegPoint, const uiint& numOfShape);
    void Calculate_J_invJ(const vvvdouble& dNdr, vector<CNode*>& vLocalNode);
    double& J(const uiint& igauss, const uiint& row, const uiint& col) {
        return mvJ[igauss][row][col];
    }
    double& inverse_J(const uiint& igauss, const uiint& row, const uiint& col) {
        return mvInvJ[igauss][row][col];
    }
    double& detJ(const uiint& igauss) {
        return mv_detJ[igauss];
    }
    double& detJ33(const vvdouble& Jmat);
    vvdouble& inverse33(const vvdouble& Jmat);
};
#endif	/* _JACOBIAN_H */
}
