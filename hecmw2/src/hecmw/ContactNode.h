/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ContactNode.h
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
#include "CommonStd.h"
#include "Vertex.h"
#include <vector>
#include <iostream>
#include "Logger.h"

namespace pmw
{
#ifndef _CONTACTNODE_H
#define	_CONTACTNODE_H
class CContactNode:public CVertex
{
public:
    CContactNode();
    virtual ~CContactNode();
protected:
    uiint myRank;
    uiint mMeshID;
    uiint mNodeID;
    uiint mLevel;

    bool mbMesh;
    bool mbNode;

    vdouble mvDisplacement;
    vdouble mvScalar;
    uiint mMFaceID;

    bool mbSlave;
    bool mbMarkingMFace;
    map<uiint, uiint, less<uiint> > mmMasterFaceID;
    vector<bool>  mvbMasterFaceID;
    vector<uiint> mvKnotID;

    bool  mbOverlap;    //--重複点判定
    vuint mvOverlapRank;//--重複点ランク:自身ランクを含む

public:
    void markingSelfMesh();//--自身に存在するMesh
    bool isSelfMesh() {
        return mbMesh;
    }
    void markingSelfNode();//--自身に存在するNode
    bool isSelfNode() {
        return mbNode;
    }

    void  setLevel(const uiint& level) {
        mLevel= level;
    }
    uiint& getLevel() {
        return mLevel;
    }
    void pushLevelMarking() {
        mvbMasterFaceID.push_back(false);
    }

    void  setRank(const uiint& rank) {
        myRank= rank;
    }
    uiint& getRank() {
        return myRank;
    }

    void  setMeshID(const uiint& meshID) {
        mMeshID= meshID;
    }
    uiint& getMeshID() {
        return mMeshID;
    }

    void  setNodeID(const uiint& nodeID) {
        mNodeID= nodeID;
    }
    uiint& getNodeID() {
        return mNodeID;
    }

    void resizeDisp(const uiint& dof);
    void initDisp();
    void setDisp(const uiint& idof, const double& disp) {
        mvDisplacement[idof]= disp;
    }
    double& getDisp(const uiint& idof) {
        return mvDisplacement[idof];
    }
    uiint getNumOfDisp() {
        return mvDisplacement.size();
    }
    void resizeScalar(const uiint& numOfScalar);
    void initScalar();
    void setScalar(const uiint& i, const double& val) {
        mvScalar[i]=val;
    }
    double& getScalar(const uiint& i) {
        return mvScalar[i];
    }
    uiint getNumOfScalar() {
        return mvScalar.size();
    }

    void markingSlave();
    bool isSlave() {
        return mbSlave;
    }

    void setMasterFaceID(const uiint& faceID, const uiint& level);
    uiint& getMasterFaceID(const uiint& level);
    bool have_MasterFaceID(const uiint& level);

    void  resizeOctreeID(const uiint& res_size);
    void  setOctreeID(const uiint& layer, const uiint& knot_id);
    uiint& getOctreeID(const uiint& layer) {
        return mvKnotID[layer];
    }

    // オーバーラップ節点:rankが重複する節点
    void  markingOverlap() {
        mbOverlap=true;
    }
    bool& isOverlap() {
        return mbOverlap;
    }
    void  addOverlapRank(const uiint& rank);//-- 重複点ランク保存:自身ランクを含む.
    uiint getOverlapMinRank();//--- AssyMatrixコンストラクタでの通信ランク決定で利用
    uiint getNumOfOverlapRank() {
        return mvOverlapRank.size();
    }
    vuint& getOverlapRank() {
        return mvOverlapRank;
    }
    uiint& getOverlapRank(const uiint& irank) {
        return mvOverlapRank[irank];
    }
    void  sort_OverlapRank();

};
#endif	/* _CONTACTNODE_H */
}
