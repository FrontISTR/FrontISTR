/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Jacobian.cpp
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
#include <vector>
#include "Vertex.h"
#include "Jacobian.h"
using namespace pmw;
CJacobian::CJacobian()
{
    uiint iaxis;
    mvJ33.resize(3);
    mvInv33.resize(3);
    for(iaxis=0; iaxis< 3; iaxis++) {
        mvJ33[iaxis].resize(3);
        mvInv33[iaxis].resize(3);
    };
}
CJacobian::~CJacobian()
{
    ;
}
void CJacobian::setupRegion(const uiint& numOfIntegPoint, const uiint& numOfShape)
{
    uiint igauss,iaxis;
    mvJ.resize(numOfIntegPoint);
    mvInvJ.resize(numOfIntegPoint);
    mv_detJ.resize(numOfIntegPoint);
    for(igauss=0; igauss< numOfIntegPoint; igauss++) {
        mvJ[igauss].resize(3);
        mvInvJ[igauss].resize(3);
        for(iaxis=0; iaxis< 3; iaxis++) {
            mvJ[igauss][iaxis].resize(3);
            mvInvJ[igauss][iaxis].resize(3);
        };
    };
}
void CJacobian::clearRegion(const uiint& numOfIntegPoint, const uiint& numOfShape)
{
    uiint igauss,iaxis;
    for(igauss=0; igauss< numOfIntegPoint; igauss++) {
        for(iaxis=0; iaxis< 3; iaxis++) {
            mvJ[igauss][iaxis].clear();
            mvInvJ[igauss][iaxis].clear();
        };
        mvJ[igauss].clear();
        mvInvJ[igauss].clear();
    };
    mvJ.clear();
    mvInvJ.clear();
    mv_detJ.clear();
}
void CJacobian::Calculate_J_invJ(const vvvdouble& dNdr, vector<CNode*>& vElemNode)
{
    uiint ig,numOfIntg;
    numOfIntg= dNdr.size();
    uiint ish,numOfShape;
    numOfShape= vElemNode.size();
    CNode* pNode;
    double X,Y,Z;
    for(ig=0; ig< numOfIntg; ig++) {
        uiint ir;
        double sum_a, sum_b, sum_c;
        for(ir=0; ir< 3; ir++) {
            sum_a=0.0;
            sum_b=0.0;
            sum_c=0.0;
            for(ish=0; ish< numOfShape; ish++) {
                pNode= vElemNode[ish];
                X= vElemNode[ish]->getX();
                Y= vElemNode[ish]->getY();
                Z= vElemNode[ish]->getZ();
                sum_a += dNdr[ig][ish][ir] * X;
                sum_b += dNdr[ig][ish][ir] * Y;
                sum_c += dNdr[ig][ish][ir] * Z;
            };
            mvJ33[0][ir]= sum_a;
            mvJ33[1][ir]= sum_b;
            mvJ33[2][ir]= sum_c;
        };
        uiint row,col;
        for(row=0; row< 3; row++) {
            for(col=0; col< 3; col++) {
                mvJ[ig][row][col]= mvJ33[row][col];
            };
        };
        mv_detJ[ig]= detJ33(mvJ33);
        inverse33(mvJ33);
        for(row=0; row< 3; row++) {
            for(col=0; col< 3; col++) {
                mvInvJ[ig][row][col]= mvInv33[row][col];
            };
        };
    };
}
double& CJacobian::detJ33(const vvdouble& Jmat)
{
    Utility::CLogger* pLogger= Utility::CLogger::Instance();
    m_detJ=  Jmat[0][0]*Jmat[1][1]*Jmat[2][2]
             +Jmat[1][0]*Jmat[2][1]*Jmat[0][2]
             +Jmat[2][0]*Jmat[0][1]*Jmat[1][2]
             -Jmat[2][0]*Jmat[1][1]*Jmat[0][2]
             -Jmat[1][0]*Jmat[0][1]*Jmat[2][2]
             -Jmat[0][0]*Jmat[2][1]*Jmat[1][2];
    if(m_detJ==0.0) pLogger->Info(Utility::LoggerMode::Error,"Math error in CJacobain, det|J|==0.0");
    return m_detJ;
}
vvdouble& CJacobian::inverse33(const vvdouble& J)
{
    double Recip= 1.0 / detJ33(J);
    mvInv33[0][0]= Recip*( J[1][1]*J[2][2] -J[2][1]*J[1][2] );
    mvInv33[0][1]= Recip*(-J[0][1]*J[2][2] +J[2][1]*J[0][2] );
    mvInv33[0][2]= Recip*( J[0][1]*J[1][2] -J[1][1]*J[0][2] );
    mvInv33[1][0]= Recip*(-J[1][0]*J[2][2] +J[2][0]*J[1][2] );
    mvInv33[1][1]= Recip*( J[0][0]*J[2][2] -J[2][0]*J[0][2] );
    mvInv33[1][2]= Recip*(-J[0][0]*J[1][2] +J[1][0]*J[0][2] );
    mvInv33[2][0]= Recip*( J[1][0]*J[2][1] -J[2][0]*J[1][1] );
    mvInv33[2][1]= Recip*(-J[0][0]*J[2][1] +J[2][0]*J[0][1] );
    mvInv33[2][2]= Recip*( J[0][0]*J[1][1] -J[1][0]*J[0][1] );
    return mvInv33;
}
