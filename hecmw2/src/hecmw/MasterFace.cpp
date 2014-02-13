/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/MasterFace.cpp
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
#include "ContactNode.h"
#include <vector>
#include "Vertex.h"
#include "MasterFace.h"
using namespace pmw;
CMasterFace::CMasterFace()
{
    mvVector.resize(3);
    mS=0.0;
    mT=0.0;
    mParam_A=0.0;
    mParam_B=0.0;
    mParam_C=0.0;
    mParam_D=0.0;
    mParam_E=0.0;
    mParam_F=0.0;
    mvParam_R.resize(3);
    mvParam_R[0]=0.0;
    mvParam_R[1]=0.0;
    mvParam_R[2]=0.0;
}
CMasterFace::~CMasterFace()
{
}
CSkinFace* CMasterFace::generateFace()
{
    mpOtherFace = new CMasterFace;
    return mpOtherFace;
}
vdouble& CMasterFace::CrossProduct(CContactNode *pConNode, const uiint& ivert, const uiint& jvert)
{
    CContactNode *pVert0,*pVert1;
    pVert0= mvConNode[ivert];
    pVert1= mvConNode[jvert];
    double X1= pVert0->getX() - pConNode->getX();
    double Y1= pVert0->getY() - pConNode->getY();
    double Z1= pVert0->getZ() - pConNode->getZ();
    double X2= pVert1->getX() - pConNode->getX();
    double Y2= pVert1->getY() - pConNode->getY();
    double Z2= pVert1->getZ() - pConNode->getZ();
    mvVector[0]= Y1*Z2 - Z1*Y2;
    mvVector[1]= Z1*X2 - X1*Z2;
    mvVector[2]= X1*Y2 - Y1*X2;
    return mvVector;
}
double CMasterFace::VectorLength(const vdouble& vec)
{
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}
double CMasterFace::VectorLength(CContactNode* pConNode0, CContactNode* pConNode1)
{
    double dX= pConNode1->getX() - pConNode0->getX();
    double dY= pConNode1->getY() - pConNode0->getY();
    double dZ= pConNode1->getZ() - pConNode0->getZ();
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(CContactNode* pConNode, const vdouble& vec)
{
    double dX= vec[0]-pConNode->getX();
    double dY= vec[1]-pConNode->getY();
    double dZ= vec[2]-pConNode->getZ();
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(const vdouble& vBasePoint, const vdouble& vPoint)
{
    double dX= vPoint[0]-vBasePoint[0];
    double dY= vPoint[1]-vBasePoint[1];
    double dZ= vPoint[2]-vBasePoint[2];
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::VectorLength(const vdouble& vec, CContactNode* pConNode)
{
    double dX= pConNode->getX() - vec[0];
    double dY= pConNode->getY() - vec[1];
    double dZ= pConNode->getZ() - vec[2];
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}
double CMasterFace::LinearInter(const double& coeff, const uiint& idof,
                                CContactNode* pConBaseNode, CContactNode* pConNode, const uiint& valType)
{
    double Value;
    switch(valType) {
    case(MPCValueType::Displacement):
        Value= coeff*(pConNode->getDisp(idof) - pConBaseNode->getDisp(idof)) + pConBaseNode->getDisp(idof);
        break;
    case(MPCValueType::Scalar):
        Value= coeff*(pConNode->getScalar(idof) - pConBaseNode->getScalar(idof)) + pConBaseNode->getScalar(idof);
        break;
    }
    return Value;
}
double CMasterFace::LinearInter(const double& coeff, const uiint& idof, const vdouble& vBaseVal, const vdouble& vVal)
{
    return coeff*(vVal[idof]-vBaseVal[idof]) + vBaseVal[idof];
}
double& CMasterFace::DotProduct()
{
    mDotValue= mvNormalVector[0]*mvVector[0] + mvNormalVector[1]*mvVector[1] + mvNormalVector[2]*mvVector[2];
    return mDotValue;
}
void CMasterFace::addSlaveNode(CContactNode* pConNode)
{
    int iprodp,iprodm;
    int nNumOfVert;
    uiint ivert,jvert;
    nNumOfVert= this->getNumOfVert();
    iprodp=0;
    iprodm=0;
    for(ivert=0; ivert< nNumOfVert; ivert++) {
        if(ivert==nNumOfVert-1) {
            jvert=0;
        } else {
            jvert= ivert+1;
        }
        CrossProduct(pConNode, ivert, jvert);
        if(DotProduct() > 0) iprodp++;
        if(DotProduct() < 0) iprodm--;
        if(DotProduct()== 0) {
            iprodp++;
            iprodm--;
        }
    };
    if(iprodp==nNumOfVert || iprodm== -nNumOfVert) {
        mvSlaveConNode.push_back(pConNode);
        mmSlaveID2Index[pConNode->getID()]= mvSlaveConNode.size()-1;
        pConNode->markingSlave();
        pConNode->setMasterFaceID(mID, mLevel);
        vdouble vCoef;
        uiint nNumOfFaceNode = mvConNode.size();
        vCoef.resize(nNumOfFaceNode);
        mvvCoef.push_back(vCoef);
    }
}
double& CMasterFace::s_NearCrossP()
{
    mS= ((mParam_B * mParam_F)-(mParam_C * mParam_E))/mParam_D;
    return mS;
}
double& CMasterFace::t_NearCrossP()
{
    mT= ((mParam_A * mParam_F)-(mParam_B * mParam_C))/mParam_D;
    return mT;
}
bool CMasterFace::NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, CContactNode* pL2P0)
{
    bool bCrossP(true);
    r_NearCrossP(L1P0, pL2P0);
    a_NearCrossP(vL1);
    b_NearCrossP(vL1,vL2);
    c_NearCrossP(vL1);
    d_NearCrossP();
    e_NearCrossP(vL2);
    f_NearCrossP(vL2);
    if( mParam_D <= 1.0E-6) bCrossP=false;
    return bCrossP;
}
bool CMasterFace::NearCrossP(const vdouble& vL1, const vdouble& vL2, CContactNode* pL1P0, CContactNode* pL2P0)
{
    bool bCrossP(true);
    r_NearCrossP(pL1P0, pL2P0);
    a_NearCrossP(vL1);
    b_NearCrossP(vL1,vL2);
    c_NearCrossP(vL1);
    d_NearCrossP();
    e_NearCrossP(vL2);
    f_NearCrossP(vL2);
    if( mParam_D <= 1.0E-6) bCrossP=false;
    return bCrossP;
}
double& CMasterFace::a_NearCrossP(const vdouble& vL1)
{
    return mParam_A= vL1[0]*vL1[0] + vL1[1]*vL1[1] + vL1[2]*vL1[2];
}
double& CMasterFace::b_NearCrossP(const vdouble& vL1, const vdouble& vL2)
{
    return mParam_B= vL1[0]*vL2[0] + vL1[1]*vL2[1] + vL1[2]*vL2[2];
}
double& CMasterFace::c_NearCrossP(const vdouble& vL1)
{
    return mParam_C= vL1[0]*mvParam_R[0] + vL1[1]*mvParam_R[1] + vL1[2]*mvParam_R[2];
}
double& CMasterFace::d_NearCrossP()
{
    return mParam_D= (mParam_A * mParam_E)-(mParam_B * mParam_B);
}
double& CMasterFace::e_NearCrossP(const vdouble& vL2)
{
    return mParam_E= vL2[0]*vL2[0] + vL2[1]*vL2[1] + vL2[2]*vL2[2];
}
double& CMasterFace::f_NearCrossP(const vdouble& vL2)
{
    return mParam_F= vL2[0]*mvParam_R[0] + vL2[1]*mvParam_R[1] + vL2[2]*mvParam_R[2];
}
vdouble& CMasterFace::r_NearCrossP(const vdouble& L1P0, const vdouble& L2P0)
{
    mvParam_R[0]= L1P0[0]-L2P0[0];
    mvParam_R[1]= L1P0[1]-L2P0[1];
    mvParam_R[2]= L1P0[2]-L2P0[2];
    return mvParam_R;
}
vdouble& CMasterFace::r_NearCrossP(CContactNode* pL1P0, CContactNode* pL2P0)
{
    mvParam_R[0]= pL1P0->getX() - pL2P0->getX();
    mvParam_R[1]= pL1P0->getY() - pL2P0->getY();
    mvParam_R[2]= pL1P0->getZ() - pL2P0->getZ();
    return mvParam_R;
}
vdouble& CMasterFace::r_NearCrossP(const vdouble& L1P0, CContactNode* pL2P0)
{
    mvParam_R[0]= L1P0[0] - pL2P0->getX();
    mvParam_R[1]= L1P0[1] - pL2P0->getY();
    mvParam_R[2]= L1P0[2] - pL2P0->getZ();
    return mvParam_R;
}
double& CMasterFace::s_ProjecP(const vdouble& vV, CContactNode* pP, CContactNode* pP3)
{
    double dX,dY,dZ;
    dX= pP->getX() - pP3->getX();
    dY= pP->getY() - pP3->getY();
    dZ= pP->getZ() - pP3->getZ();
    mT= vV[0]*dX + vV[1]*dY + vV[2]*dZ;
    return mS;
}
double& CMasterFace::t_ProjecP(const vdouble& vV, CContactNode* pP, CContactNode* pP2)
{
    double dX,dY,dZ;
    dX= pP->getX() - pP2->getX();
    dY= pP->getY() - pP2->getY();
    dZ= pP->getZ() - pP2->getZ();
    mT= vV[0]*dX + vV[1]*dY + vV[2]*dZ;
    return mT;
}
void CMasterFace::Normalized(vdouble& vV)
{
    double dLength= sqrt(vV[0]*vV[0] + vV[1]*vV[1] + vV[2]*vV[2]);
    vV[0]/= dLength;
    vV[1]/= dLength;
    vV[2]/= dLength;
}
double CMasterFace::CoefTermA(CContactNode* pOpposNode, const vdouble& inP, CContactNode* pEdgeNode0, CContactNode* pEdgeNode1)
{
    double nume;
    double deno;
    double termA;
    nume= VectorLength(pOpposNode, inP);
    deno= VectorLength(pEdgeNode0,pEdgeNode1);
    termA= nume/deno;
    return termA;
}
double CMasterFace::CoefTermB(const vdouble& inP, CContactNode* pSlaveP, const vdouble& inP0, const vdouble& inP1)
{
    double nume;
    double deno;
    double termB;
    nume= VectorLength(pSlaveP, inP);
    deno= VectorLength(inP0, inP1);
    termB= nume/deno;
    return termB;
}
void CMasterFace::CalcSlave(const uiint& islave, const uiint& valType)
{
    CContactNode* pSlaveP= mvSlaveConNode[islave];
    uiint numOfVert= this->getNumOfVert();
    vdouble vOP;
    vOP.resize(3);
    vdouble vNP0,vNP1;
    vNP0.resize(3);
    vNP1.resize(3);
    vdouble vE,vF;
    vE.resize(3);
    vF.resize(3);
    vdouble vPOP;
    vPOP.resize(3);
    vdouble vP2P3;
    vP2P3.resize(3);
    vdouble vP1P0;
    vP1P0.resize(3);
    vdouble vInterP0,vInterP1;
    vInterP0.resize(3);
    vInterP1.resize(3);
    double dLenP2P3,dLenP2InP,dLenP0P1,dLenP1InP;
    vdouble vValInterP0,vValInterP1;
    uiint numOfDOF,idof;
    double dLenInP0InP1,dLenInP0SlaveP;
    vdouble vValSlaveP;
    double  aArea;
    vdouble vArea;
    vArea.resize(3);
    switch(valType) {
    case(MPCValueType::Displacement):
        numOfDOF=3;
        vValInterP0.resize(numOfDOF);
        vValInterP1.resize(numOfDOF);
        vValSlaveP.resize(numOfDOF);
        break;
    case(MPCValueType::Scalar):
        numOfDOF=mvConNode[0]->getNumOfScalar();
        vValInterP0.resize(numOfDOF);
        vValInterP1.resize(numOfDOF);
        vValSlaveP.resize(numOfDOF);
        break;
    }
    bool bEF;
    uiint ideb;
    double deb;
    switch(numOfVert) {
    case(4):
        vE[0]= mvConNode[0]->getX() - mvConNode[3]->getX();
        vE[1]= mvConNode[0]->getY() - mvConNode[3]->getY();
        vE[2]= mvConNode[0]->getZ() - mvConNode[3]->getZ();
        vF[0]= mvConNode[1]->getX() - mvConNode[2]->getX();
        vF[1]= mvConNode[1]->getY() - mvConNode[2]->getY();
        vF[2]= mvConNode[1]->getZ() - mvConNode[2]->getZ();
        Normalized(vE);
        Normalized(vF);
        bEF= NearCrossP(vE,vF, mvConNode[3],mvConNode[2]);
        if(bEF) {
            s_NearCrossP();
            t_NearCrossP();
            vNP0[0]= mS*vE[0] + mvConNode[3]->getX();
            vNP0[1]= mS*vE[1] + mvConNode[3]->getY();
            vNP0[2]= mS*vE[2] + mvConNode[3]->getZ();
            vNP1[0]= mT*vF[0] + mvConNode[2]->getX();
            vNP1[1]= mT*vF[1] + mvConNode[2]->getY();
            vNP1[2]= mT*vF[2] + mvConNode[2]->getZ();
            vOP[0]= (vNP0[0]+vNP1[0])*0.5;
            vOP[1]= (vNP0[1]+vNP1[1])*0.5;
            vOP[2]= (vNP0[2]+vNP1[2])*0.5;
            vPOP[0]= pSlaveP->getX() - vOP[0];
            vPOP[1]= pSlaveP->getY() - vOP[1];
            vPOP[2]= pSlaveP->getZ() - vOP[2];
        } else {
            vPOP[0]= vE[0];
            vPOP[1]= vE[1];
            vPOP[2]= vE[2];
            vOP[0]= pSlaveP->getX();
            vOP[1]= pSlaveP->getY();
            vOP[2]= pSlaveP->getZ();
        }
        vP2P3[0]= mvConNode[3]->getX() - mvConNode[2]->getX();
        vP2P3[1]= mvConNode[3]->getY() - mvConNode[2]->getY();
        vP2P3[2]= mvConNode[3]->getZ() - mvConNode[2]->getZ();
        vP1P0[0]= mvConNode[0]->getX() - mvConNode[1]->getX();
        vP1P0[1]= mvConNode[0]->getY() - mvConNode[1]->getY();
        vP1P0[2]= mvConNode[0]->getZ() - mvConNode[1]->getZ();
        Normalized(vPOP);
        Normalized(vP2P3);
        Normalized(vP1P0);
        NearCrossP(vPOP,vP2P3, vOP,mvConNode[2]);
        s_NearCrossP();
        t_NearCrossP();
        vNP0[0]= mS*vPOP[0] + vOP[0];
        vNP0[1]= mS*vPOP[1] + vOP[1];
        vNP0[2]= mS*vPOP[2] + vOP[2];
        vNP1[0]= mT*vP2P3[0] + mvConNode[2]->getX();
        vNP1[1]= mT*vP2P3[1] + mvConNode[2]->getY();
        vNP1[2]= mT*vP2P3[2] + mvConNode[2]->getZ();
        vInterP0[0]= (vNP0[0]+vNP1[0])*0.5;
        vInterP0[1]= (vNP0[1]+vNP1[1])*0.5;
        vInterP0[2]= (vNP0[2]+vNP1[2])*0.5;
        dLenP2P3= VectorLength(mvConNode[2],mvConNode[3]);
        dLenP2InP= VectorLength(mvConNode[2],vInterP0);
        for(idof=0; idof < numOfDOF; idof++) {
            vValInterP0[idof]= LinearInter(dLenP2InP/dLenP2P3,idof,mvConNode[2],mvConNode[3],valType);
        };
        NearCrossP(vPOP,vP1P0, vOP,mvConNode[1]);
        s_NearCrossP();
        t_NearCrossP();
        vNP0[0]= mS*vPOP[0] + vOP[0];
        vNP0[1]= mS*vPOP[1] + vOP[1];
        vNP0[2]= mS*vPOP[2] + vOP[2];
        vNP1[0]= mT*vP1P0[0] + mvConNode[1]->getX();
        vNP1[1]= mT*vP1P0[1] + mvConNode[1]->getY();
        vNP1[2]= mT*vP1P0[2] + mvConNode[1]->getZ();
        vInterP1[0]= (vNP0[0]+vNP1[0])*0.5;
        vInterP1[1]= (vNP0[1]+vNP1[1])*0.5;
        vInterP1[2]= (vNP0[2]+vNP1[2])*0.5;
        dLenP0P1= VectorLength(mvConNode[1],mvConNode[0]);
        dLenP1InP= VectorLength(mvConNode[1],vInterP1);
        for(idof=0; idof < numOfDOF; idof++) {
            vValInterP1[idof]= LinearInter(dLenP1InP/dLenP0P1,idof,mvConNode[1],mvConNode[0],valType);
        };
        mvvCoef[islave][0]= CoefTermA(mvConNode[1],vInterP1,mvConNode[0],mvConNode[1]) * CoefTermB(vInterP0,pSlaveP, vInterP1,vInterP0);
        mvvCoef[islave][1]= CoefTermA(mvConNode[0],vInterP1,mvConNode[0],mvConNode[1]) * CoefTermB(vInterP0,pSlaveP, vInterP1,vInterP0);
        mvvCoef[islave][2]= CoefTermA(mvConNode[3],vInterP0,mvConNode[2],mvConNode[3]) * CoefTermB(vInterP1,pSlaveP, vInterP1,vInterP0);
        mvvCoef[islave][3]= CoefTermA(mvConNode[2],vInterP0,mvConNode[2],mvConNode[3]) * CoefTermB(vInterP1,pSlaveP, vInterP1,vInterP0);
        if(mnOrder==ElementOrder::Second) {
            mvvCoef[islave][4] = 0.0;
            mvvCoef[islave][5] = 0.0;
            mvvCoef[islave][6] = 0.0;
            mvvCoef[islave][7] = 0.0;
        }
        break;
    case(3):
        vOP[0]= mvConNode[2]->getX();
        vOP[1]= mvConNode[2]->getY();
        vOP[2]= mvConNode[2]->getZ();
        vPOP[0]= pSlaveP->getX() - vOP[0];
        vPOP[1]= pSlaveP->getY() - vOP[1];
        vPOP[2]= pSlaveP->getZ() - vOP[2];
        vP1P0[0]= mvConNode[0]->getX() - mvConNode[1]->getX();
        vP1P0[1]= mvConNode[0]->getY() - mvConNode[1]->getY();
        vP1P0[2]= mvConNode[0]->getZ() - mvConNode[1]->getZ();
        Normalized(vPOP);
        Normalized(vP1P0);
        NearCrossP(vPOP,vP1P0, vOP,mvConNode[1]);
        s_NearCrossP();
        t_NearCrossP();
        vNP0[0]= mS*vPOP[0] + vOP[0];
        vNP0[1]= mS*vPOP[1] + vOP[1];
        vNP0[2]= mS*vPOP[2] + vOP[2];
        vNP1[0]= mT*vP1P0[0] + mvConNode[1]->getX();
        vNP1[1]= mT*vP1P0[1] + mvConNode[1]->getY();
        vNP1[2]= mT*vP1P0[2] + mvConNode[1]->getZ();
        vInterP1[0]= (vNP0[0]+vNP1[0])*0.5;
        vInterP1[1]= (vNP0[1]+vNP1[1])*0.5;
        vInterP1[2]= (vNP0[2]+vNP1[2])*0.5;
        dLenP0P1= VectorLength(mvConNode[1],mvConNode[0]);
        dLenP1InP= VectorLength(mvConNode[1],vInterP1);
        for(idof=0; idof < numOfDOF; idof++) {
            vValInterP1[idof]= LinearInter(dLenP1InP/dLenP0P1,idof,mvConNode[1],mvConNode[0],valType);
        };
        vInterP0[0]= mvConNode[2]->getX();
        vInterP0[1]= mvConNode[2]->getY();
        vInterP0[2]= mvConNode[2]->getZ();
        switch(valType) {
        case(MPCValueType::Displacement):
            for(idof=0; idof < numOfDOF; idof++) vValInterP0[idof]= mvConNode[2]->getDisp(idof);
            break;
        case(MPCValueType::Scalar):
            for(idof=0; idof < numOfDOF; idof++) vValInterP0[idof]= mvConNode[2]->getScalar(idof);
            break;
        }
        CrossProduct(mvConNode[0], 1,2);
        aArea= VectorLength(mvVector)*0.5;
        CrossProduct(pSlaveP, 1,2);
        vArea[0]= VectorLength(mvVector)*0.5;
        CrossProduct(pSlaveP, 2,0);
        vArea[1]= VectorLength(mvVector)*0.5;
        CrossProduct(pSlaveP, 0,1);
        vArea[2]= VectorLength(mvVector)*0.5;
        mvvCoef[islave][0]= vArea[0]/aArea;
        mvvCoef[islave][1]= vArea[1]/aArea;
        mvvCoef[islave][2]= vArea[2]/aArea;
        if(mnOrder==ElementOrder::Second) {
            mvvCoef[islave][3] = 0.0;
            mvvCoef[islave][4] = 0.0;
            mvvCoef[islave][5] = 0.0;
        }
        break;
    case(2):
        vInterP0[0]= mvConNode[0]->getX();
        vInterP0[1]= mvConNode[0]->getY();
        vInterP0[2]= mvConNode[0]->getZ();
        vInterP1[0]= mvConNode[1]->getX();
        vInterP1[1]= mvConNode[1]->getY();
        vInterP1[2]= mvConNode[1]->getZ();
        switch(valType) {
        case(MPCValueType::Displacement):
            for(idof=0; idof < numOfDOF; idof++) {
                vValInterP0[idof]= mvConNode[0]->getDisp(idof);
                vValInterP1[idof]= mvConNode[1]->getDisp(idof);
            };
            break;
        case(MPCValueType::Scalar):
            for(idof=0; idof < numOfDOF; idof++) {
                vValInterP0[idof]= mvConNode[0]->getScalar(idof);
                vValInterP1[idof]= mvConNode[1]->getScalar(idof);
            };
            break;
        }
        mvvCoef[islave][0]= CoefTermB(vInterP1,pSlaveP, vInterP0,vInterP1);
        mvvCoef[islave][1]= CoefTermB(vInterP0,pSlaveP, vInterP0,vInterP1);
        if(mnOrder==ElementOrder::Second) {
            mvvCoef[islave][2] = 0.0;
        }
        break;
    default:
        break;
    }
    dLenInP0InP1= VectorLength(vInterP0, vInterP1);
    dLenInP0SlaveP= VectorLength(vInterP0, pSlaveP);
    for(idof=0; idof < numOfDOF; idof++) {
        vValSlaveP[idof]= LinearInter(dLenInP0SlaveP/dLenInP0InP1,idof, vValInterP0,vValInterP1);
    };
    switch(valType) {
    case(MPCValueType::Displacement):
        for(idof=0; idof < numOfDOF; idof++)
            pSlaveP->setDisp(idof, vValSlaveP[idof]);
        break;
    case(MPCValueType::Scalar):
        for(idof=0; idof < numOfDOF; idof++)
            pSlaveP->setScalar(idof, vValSlaveP[idof]);
        break;
    }
}
double& CMasterFace::getCoef(const uiint& slaveID, const uiint& ivert)
{
    uiint islave= mmSlaveID2Index[slaveID];
    return mvvCoef[islave][ivert];
}
