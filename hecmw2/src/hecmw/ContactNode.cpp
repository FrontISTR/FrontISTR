/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/ContactNode.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "ContactNode.h"
using namespace pmw;
CContactNode::CContactNode()
{
    mbMesh=false;
    mbNode=false;
    mbSlave=false;
    mbMarkingMFace=false;
}
CContactNode::~CContactNode()
{
    ;
}
void CContactNode::markingSelfMesh()
{
    mbMesh=true;
}
void CContactNode::markingSelfNode()
{
    mbNode=true;
}
void CContactNode::resizeDisp(const uiint& dof)
{
    mvDisplacement.resize(dof);
}
void CContactNode::initDisp()
{
    uiint numOfDOF;
    uiint idof;
    numOfDOF= mvDisplacement.size();
    for(idof=0; idof< numOfDOF; idof++){
        mvDisplacement[idof]=0.0;
    };
}
void CContactNode::resizeScalar(const uiint& numOfScalar)
{
    mvScalar.resize(numOfScalar);
}
void CContactNode::initScalar()
{
    uiint numOfScalar= mvScalar.size();
    uiint i;
    for(i=0; i< numOfScalar; i++){
        mvScalar[i]= 0.0;
    };
}
void CContactNode::markingSlave()
{
    mbSlave=true;
}
void CContactNode::setMasterFaceID(const uiint& faceID, const uiint& level)
{
    mmMasterFaceID[level]= faceID;
    mvbMasterFaceID[level-mLevel]=true;
}
uiint& CContactNode::getMasterFaceID(const uiint& level)
{
    return mmMasterFaceID[level];
}
bool CContactNode::have_MasterFaceID(const uiint& level)
{
    return mvbMasterFaceID[level - mLevel];
}
void CContactNode::resizeOctreeID(const uiint& res_size)
{
    mvKnotID.resize(res_size);
}
void CContactNode::setOctreeID(const uiint& layer, const uiint& knot_id)
{
    mvKnotID[layer]= knot_id;
}
