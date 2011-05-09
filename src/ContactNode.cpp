/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ContactNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
void CContactNode::resizeDisp(const uint& dof)
{
    mvDisplacement.resize(dof);
}
void CContactNode::initDisp()
{
    uint numOfDOF;
    uint idof;
    numOfDOF= mvDisplacement.size();
    for(idof=0; idof< numOfDOF; idof++){
        mvDisplacement[idof]=0.0;
    };
}
void CContactNode::resizeScalar(const uint& numOfScalar)
{
    mvScalar.resize(numOfScalar);
}
void CContactNode::initScalar()
{
    uint numOfScalar= mvScalar.size();
    uint i;
    for(i=0; i< numOfScalar; i++){
        mvScalar[i]= 0.0;
    };
}
void CContactNode::markingSlave()
{
    mbSlave=true;
}
void CContactNode::setMasterFaceID(const uint& faceID, const uint& level)
{
    mmMasterFaceID[level]= faceID;
    mvbMasterFaceID[level-mLevel]=true;
}
uint& CContactNode::getMasterFaceID(const uint& level)
{
    return mmMasterFaceID[level];
}
bool CContactNode::have_MasterFaceID(const uint& level)
{
    return mvbMasterFaceID[level - mLevel];
}
void CContactNode::resizeOctreeID(const uint& res_size)
{
    mvKnotID.resize(res_size);
}
void CContactNode::setOctreeID(const uint& layer, const uint& knot_id)
{
    mvKnotID[layer]= knot_id;
}
