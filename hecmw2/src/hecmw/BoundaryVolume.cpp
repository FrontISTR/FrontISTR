/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BoundaryVolume.cpp
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
#include "EdgeTree.h"
#include "BoundaryParts.h"
#include "Element.h"
#include "BoundaryNode.h"
#include <vector>
#include "Logger.h"
#include "BoundaryVolume.h"
#include "ElementType.h"
#include "ElementProperty.h"
using namespace pmw;
CBoundaryVolume::CBoundaryVolume()
{
    ;
}
CBoundaryVolume::~CBoundaryVolume()
{
    ;
}
void CBoundaryVolume::setEdgeNeibVol(const uiint& iedge, const uiint& neibVolID)
{
    mvEdgeNeibVol[iedge].push_back(neibVolID);
}
void CBoundaryVolume::setFaceNeibVol(const uiint& iface, const uiint& neibVolID)
{
    mvFaceNeibVol[iface]= neibVolID;
}
void CBoundaryVolume::setEdgeBNode(const uiint& iedge, CBoundaryNode* pBNode)
{
    mvEdgeBNode[iedge]= pBNode;
}
void CBoundaryVolume::setFaceBNode(const uiint& iface, CBoundaryNode* pBNode)
{
    mvFaceBNode[iface]= pBNode;
}
void CBoundaryVolume::setVolBNode(CBoundaryNode* pBNode)
{
    mpVolBNode= pBNode;
}
void CBoundaryVolume::distValue(CBoundaryVolume* pProgVol, const double& coef, const vuint& vDOF)
{
    uiint idof, dof;
    double progValue;
    for(idof=0; idof < vDOF.size(); idof++) {
        dof = vDOF[idof];
        progValue = mmValue[dof]*coef;
        pProgVol->setBndValue(dof,progValue);
    };
}
uiint CBoundaryVolume::dividHexa(const uiint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog) {
    case(0):
        pProgVol->setBNode(0, mvBNode[0]);
        pProgVol->setBNode(1, mvEdgeBNode[0]);
        pProgVol->setBNode(2, mvFaceBNode[0]);
        pProgVol->setBNode(3, mvEdgeBNode[3]);
        pProgVol->setBNode(4, mvEdgeBNode[8]);
        pProgVol->setBNode(5, mvFaceBNode[4]);
        pProgVol->setBNode(6, mpVolBNode);
        pProgVol->setBNode(7, mvFaceBNode[3]);
        return 0;
        break;
    case(1):
        pProgVol->setBNode(0, mvEdgeBNode[0]);
        pProgVol->setBNode(1, mvBNode[1]);
        pProgVol->setBNode(2, mvEdgeBNode[1]);
        pProgVol->setBNode(3, mvFaceBNode[0]);
        pProgVol->setBNode(4, mvFaceBNode[4]);
        pProgVol->setBNode(5, mvEdgeBNode[9]);
        pProgVol->setBNode(6, mvFaceBNode[2]);
        pProgVol->setBNode(7, mpVolBNode);
        return 1;
        break;
    case(2):
        pProgVol->setBNode(0, mvEdgeBNode[8]);
        pProgVol->setBNode(1, mvFaceBNode[4]);
        pProgVol->setBNode(2, mpVolBNode);
        pProgVol->setBNode(3, mvFaceBNode[3]);
        pProgVol->setBNode(4, mvBNode[4]);
        pProgVol->setBNode(5, mvEdgeBNode[4]);
        pProgVol->setBNode(6, mvFaceBNode[1]);
        pProgVol->setBNode(7, mvEdgeBNode[7]);
        return 4;
        break;
    case(3):
        pProgVol->setBNode(0, mvFaceBNode[4]);
        pProgVol->setBNode(1, mvEdgeBNode[9]);
        pProgVol->setBNode(2, mvFaceBNode[2]);
        pProgVol->setBNode(3, mpVolBNode);
        pProgVol->setBNode(4, mvEdgeBNode[4]);
        pProgVol->setBNode(5, mvBNode[5]);
        pProgVol->setBNode(6, mvEdgeBNode[5]);
        pProgVol->setBNode(7, mvFaceBNode[1]);
        return 5;
        break;
    case(4):
        pProgVol->setBNode(0, mvEdgeBNode[3]);
        pProgVol->setBNode(1, mvFaceBNode[0]);
        pProgVol->setBNode(2, mvEdgeBNode[2]);
        pProgVol->setBNode(3, mvBNode[3]);
        pProgVol->setBNode(4, mvFaceBNode[3]);
        pProgVol->setBNode(5, mpVolBNode);
        pProgVol->setBNode(6, mvFaceBNode[5]);
        pProgVol->setBNode(7, mvEdgeBNode[11]);
        return 3;
        break;
    case(5):
        pProgVol->setBNode(0, mvFaceBNode[0]);
        pProgVol->setBNode(1, mvEdgeBNode[1]);
        pProgVol->setBNode(2, mvBNode[2]);
        pProgVol->setBNode(3, mvEdgeBNode[2]);
        pProgVol->setBNode(4, mpVolBNode);
        pProgVol->setBNode(5, mvFaceBNode[2]);
        pProgVol->setBNode(6, mvEdgeBNode[10]);
        pProgVol->setBNode(7, mvFaceBNode[5]);
        return 2;
        break;
    case(6):
        pProgVol->setBNode(0, mvFaceBNode[3]);
        pProgVol->setBNode(1, mpVolBNode);
        pProgVol->setBNode(2, mvFaceBNode[5]);
        pProgVol->setBNode(3, mvEdgeBNode[11]);
        pProgVol->setBNode(4, mvEdgeBNode[7]);
        pProgVol->setBNode(5, mvFaceBNode[1]);
        pProgVol->setBNode(6, mvEdgeBNode[6]);
        pProgVol->setBNode(7, mvBNode[7]);
        return 7;
        break;
    case(7):
        pProgVol->setBNode(0, mpVolBNode);
        pProgVol->setBNode(1, mvFaceBNode[2]);
        pProgVol->setBNode(2, mvEdgeBNode[10]);
        pProgVol->setBNode(3, mvFaceBNode[5]);
        pProgVol->setBNode(4, mvFaceBNode[1]);
        pProgVol->setBNode(5, mvEdgeBNode[5]);
        pProgVol->setBNode(6, mvBNode[6]);
        pProgVol->setBNode(7, mvEdgeBNode[6]);
        return 6;
        break;
    }
}
uiint CBoundaryVolume::dividTetra(const uiint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog) {
    case(0):
        pProgVol->setBNode(0, mvEdgeBNode[2]);
        pProgVol->setBNode(1, mvBNode[0]);
        pProgVol->setBNode(2, mvEdgeBNode[0]);
        pProgVol->setBNode(3, mvFaceBNode[0]);
        pProgVol->setBNode(4, mvFaceBNode[3]);
        pProgVol->setBNode(5, mvEdgeBNode[3]);
        pProgVol->setBNode(6, mvFaceBNode[1]);
        pProgVol->setBNode(7, mpVolBNode);
        return 0;
        break;
    case(1):
        pProgVol->setBNode(0, mvBNode[2]);
        pProgVol->setBNode(1, mvEdgeBNode[2]);
        pProgVol->setBNode(2, mvFaceBNode[0]);
        pProgVol->setBNode(3, mvEdgeBNode[1]);
        pProgVol->setBNode(4, mvEdgeBNode[5]);
        pProgVol->setBNode(5, mvFaceBNode[3]);
        pProgVol->setBNode(6, mpVolBNode);
        pProgVol->setBNode(7, mvFaceBNode[2]);
        return 2;
        break;
    case(2):
        pProgVol->setBNode(0, mvFaceBNode[0]);
        pProgVol->setBNode(1, mvEdgeBNode[0]);
        pProgVol->setBNode(2, mvBNode[1]);
        pProgVol->setBNode(3, mvEdgeBNode[1]);
        pProgVol->setBNode(4, mpVolBNode);
        pProgVol->setBNode(5, mvFaceBNode[1]);
        pProgVol->setBNode(6, mvEdgeBNode[4]);
        pProgVol->setBNode(7, mvFaceBNode[2]);
        return 1;
        break;
    case(3):
        pProgVol->setBNode(0, mvFaceBNode[3]);
        pProgVol->setBNode(1, mvEdgeBNode[3]);
        pProgVol->setBNode(2, mvFaceBNode[1]);
        pProgVol->setBNode(3, mpVolBNode);
        pProgVol->setBNode(4, mvEdgeBNode[5]);
        pProgVol->setBNode(5, mvBNode[3]);
        pProgVol->setBNode(6, mvEdgeBNode[4]);
        pProgVol->setBNode(7, mvFaceBNode[2]);
        return 3;
        break;
    }
}
uiint CBoundaryVolume::dividPrism(const uiint& iprog, CBoundaryVolume* pProgVol)
{
    switch(iprog) {
    case(0):
        pProgVol->setBNode(0, mvBNode[2]);
        pProgVol->setBNode(1, mvEdgeBNode[1]);
        pProgVol->setBNode(2, mvFaceBNode[0]);
        pProgVol->setBNode(3, mvEdgeBNode[2]);
        pProgVol->setBNode(4, mvEdgeBNode[5]);
        pProgVol->setBNode(5, mvFaceBNode[4]);
        pProgVol->setBNode(6, mpVolBNode);
        pProgVol->setBNode(7, mvFaceBNode[3]);
        return 2;
        break;
    case(1):
        pProgVol->setBNode(0, mvEdgeBNode[1]);
        pProgVol->setBNode(1, mvBNode[0]);
        pProgVol->setBNode(2, mvEdgeBNode[0]);
        pProgVol->setBNode(3, mvFaceBNode[0]);
        pProgVol->setBNode(4, mvFaceBNode[4]);
        pProgVol->setBNode(5, mvEdgeBNode[3]);
        pProgVol->setBNode(6, mvFaceBNode[2]);
        pProgVol->setBNode(7, mpVolBNode);
        return 0;
        break;
    case(2):
        pProgVol->setBNode(0, mvFaceBNode[0]);
        pProgVol->setBNode(1, mvEdgeBNode[0]);
        pProgVol->setBNode(2, mvBNode[1]);
        pProgVol->setBNode(3, mvEdgeBNode[2]);
        pProgVol->setBNode(4, mpVolBNode);
        pProgVol->setBNode(5, mvFaceBNode[2]);
        pProgVol->setBNode(6, mvEdgeBNode[4]);
        pProgVol->setBNode(7, mvFaceBNode[3]);
        return 1;
        break;
    case(3):
        pProgVol->setBNode(0, mvEdgeBNode[5]);
        pProgVol->setBNode(1, mvFaceBNode[4]);
        pProgVol->setBNode(2, mpVolBNode);
        pProgVol->setBNode(3, mvFaceBNode[3]);
        pProgVol->setBNode(4, mvBNode[5]);
        pProgVol->setBNode(5, mvEdgeBNode[8]);
        pProgVol->setBNode(6, mvFaceBNode[1]);
        pProgVol->setBNode(7, mvEdgeBNode[7]);
        return 5;
        break;
    case(4):
        pProgVol->setBNode(0, mvFaceBNode[4]);
        pProgVol->setBNode(1, mvEdgeBNode[3]);
        pProgVol->setBNode(2, mvFaceBNode[2]);
        pProgVol->setBNode(3, mpVolBNode);
        pProgVol->setBNode(4, mvEdgeBNode[8]);
        pProgVol->setBNode(5, mvBNode[3]);
        pProgVol->setBNode(6, mvEdgeBNode[6]);
        pProgVol->setBNode(7, mvFaceBNode[1]);
        return 3;
        break;
    case(5):
        pProgVol->setBNode(0, mpVolBNode);
        pProgVol->setBNode(1, mvFaceBNode[2]);
        pProgVol->setBNode(2, mvEdgeBNode[4]);
        pProgVol->setBNode(3, mvFaceBNode[3]);
        pProgVol->setBNode(4, mvFaceBNode[1]);
        pProgVol->setBNode(5, mvEdgeBNode[6]);
        pProgVol->setBNode(6, mvBNode[4]);
        pProgVol->setBNode(7, mvEdgeBNode[7]);
        return 4;
        break;
    }
}
double CBoundaryVolume::tetraVolume(CNode* pNode0, CNode* pNode1, CNode* pNode2, CNode* pNode3)
{
    double vol(0.0);
    double u[3],v[3],w[3];
    u[0]= pNode1->getX() - pNode0->getX();
    u[1]= pNode1->getY() - pNode0->getY();
    u[2]= pNode1->getZ() - pNode0->getZ();
    v[0]= pNode2->getX() - pNode0->getX();
    v[1]= pNode2->getY() - pNode0->getY();
    v[2]= pNode2->getZ() - pNode0->getZ();
    w[0]= pNode3->getX() - pNode0->getX();
    w[1]= pNode3->getY() - pNode0->getY();
    w[2]= pNode3->getZ() - pNode0->getZ();
    double cross[3];
    cross[0] = v[1]*w[2] - v[2]*w[1];
    cross[1] = v[2]*w[0] - v[0]*w[2];
    cross[2] = v[0]*w[1] - v[1]*w[0];
    double scalar;
    scalar = u[0]*cross[0] + u[1]*cross[1] + u[2]*cross[2];
    if(scalar < 0.0) scalar *= -1.0;
    vol = (1.0/6.0)*scalar;
    return vol;
}
