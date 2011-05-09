/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   MasterFace.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef MASTER_FACE_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC
#define MASTER_FACE_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC
#include "SkinFace.h"
#include "MPCValueType.h"
#include <map> 
namespace pmw{
class CMasterFace:public CSkinFace{
public:
    CMasterFace();
    virtual ~CMasterFace();
private:
    double mDotValue;
    vdouble mvVector;
    vdouble& CrossProduct(CContactNode* pConNode,const uint& ivert, const uint& jvert);
    double&  DotProduct();
    double VectorLength(const vdouble& vec);
    double VectorLength(CContactNode *pConNode0, CContactNode *pConNode1);
    double VectorLength(CContactNode *pConNode, const vdouble& vec);
    double VectorLength(const vdouble& vBasePoint, const vdouble& vPoint);
    double VectorLength(const vdouble& vec, CContactNode *pConNode);
    double LinearInter(const double& coeff, const uint& idof, CContactNode *pConBaseNode, CContactNode *pConNode ,const uint& valType);
    double LinearInter(const double& coeff, const uint& idof, const vdouble& vBaseVal, const vdouble& vVal);
    double mS,mT, mParam_A, mParam_B, mParam_C, mParam_D, mParam_E, mParam_F;
    vdouble mvParam_R;
    double& s_NearCrossP();
    double& t_NearCrossP();
    bool NearCrossP(const vdouble& vL1, const vdouble& vL2, CContactNode* pL1P0, CContactNode* pL2P0);
    bool NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, CContactNode* pL2P0);
    double& a_NearCrossP(const vdouble& vL1);
    double& b_NearCrossP(const vdouble& vL1, const vdouble& vL2);
    double& c_NearCrossP(const vdouble& vL1);
    double& d_NearCrossP();                  
    double& e_NearCrossP(const vdouble& vL2);
    double& f_NearCrossP(const vdouble& vL2);
    vdouble& r_NearCrossP(const vdouble& L1P0, const vdouble& L2P0);
    vdouble& r_NearCrossP(CContactNode *pL1P0, CContactNode *pL2P0);
    vdouble& r_NearCrossP(const vdouble& L1P0, CContactNode *pL2P0);
    double& s_ProjecP(const vdouble& vV, CContactNode *pP,CContactNode *pP3);
    double& t_ProjecP(const vdouble& vV, CContactNode *pP,CContactNode *pP2);
    void Normalized(vdouble& vV);
    vvdouble mvvCoef;
    double CoefTermA(CContactNode *pOpposNode, const vdouble& inP, CContactNode *pEdgeNode0, CContactNode *pEdgeNode1);
    double CoefTermB(const vdouble& inP,CContactNode *pSlaveP, const vdouble& inP0, const vdouble& inP1);
    std::map<uint,uint,less<uint> > mmSlaveID2Index;
protected:
    vector<CContactNode*> mvSlaveConNode;
    virtual CSkinFace* generateFace();
public:
    virtual const char* getName(){ return "MasterFace";}
    virtual void addSlaveNode(CContactNode* pConNode);
    virtual CContactNode* getSlaveNode(const uint& index){ return mvSlaveConNode[index];}
    virtual uint getNumOfSlaveNode(){ return mvSlaveConNode.size();}
    virtual void CalcSlave(const uint& islave, const uint& valType);
    virtual double& getCoef(const uint& slaveID, const uint& ivert);
};
}
#endif
