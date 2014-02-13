/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ContactNode.cpp
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
using namespace pmw;
CContactNode::CContactNode()
{
    mbMesh=false;
    mbNode=false;
    mbSlave=false;
    mbMarkingMFace=false;

    mbOverlap=false;
}
CContactNode::~CContactNode()
{
    ;
}
//自身に存在するMesh
void CContactNode::markingSelfMesh()
{
    mbMesh=true;
}
//自身に存在するNode
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
    for(idof=0; idof< numOfDOF; idof++) {
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
    for(i=0; i< numOfScalar; i++) {
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

//--
// オーバーラップランク点の保存:自身のランクも含む
//--
void CContactNode::addOverlapRank(const uiint& rank)
{
    uiint nNumRank= mvOverlapRank.size();

    // 最初のランク
    if(nNumRank==0) {
        mvOverlapRank.push_back(rank);
        return;
    }

    // 一個以上のランクが既にある場合
    bool bFind(false);
    for(uiint irank=0; irank < nNumRank; irank++) {
        if(mvOverlapRank[irank]==rank) bFind=true;
    };
    if(!bFind) mvOverlapRank.push_back(rank);

}
//--
// sort for オーバーラップランク
//--
void CContactNode::sort_OverlapRank()
{
    std::sort(mvOverlapRank.begin(), mvOverlapRank.end());
}
//--
// オーバーラップ点の最小Rank : AssyMatrixコンストラクターでの通信ランク決定に利用.
//--
uiint CContactNode::getOverlapMinRank()
{
    if(!mbOverlap) {
        Utility::CLogger* pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Warn, "No Overlap Node, ContactNode::getOverlapMinRank");

        return IINT_MAX;
    }

    uiint nNumOfOverlap = mvOverlapRank.size();
    uiint nMinRank;
    if(nNumOfOverlap > 0) {
        // 初期値
        nMinRank= mvOverlapRank[0];
    } else {
        // Error
        Utility::CLogger* pLogger= Utility::CLogger::Instance();
        pLogger->Info(Utility::LoggerMode::Error, "Did not set OverlapRank, ContactNode::getOverlapMinRank");

        return IINT_MAX;
    }
    //--
    // 最小値
    //--
    for(uiint i=0; i < nNumOfOverlap; i++) {
        if(nMinRank > mvOverlapRank[i]) nMinRank= mvOverlapRank[i];
    };

    return nMinRank;
}


