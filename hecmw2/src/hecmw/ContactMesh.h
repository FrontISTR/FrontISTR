/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ContactMesh.h
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
#ifndef CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F
#define CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F
#include "SkinFace.h"
#include "MasterFace.h"
#include "BoundingBox.h"
#include "MPCValueType.h"
#include "OctreeKnot.h"
#include "Film.h"

#include <map>
#include <utility>
typedef std::pair<pmw::CContactNode*, pmw::CContactNode*> PairConNode;
namespace pmw
{
class CContactMesh
{
public:
    CContactMesh();
    virtual ~CContactMesh();
protected:
    uiint mContactID;
    uiint mLevel;
    uiint myRank;

    uiint mnProp;
    CFilm *mpFilm;//--------- 熱伝達係数 : 線形方程式番号に対応

    vuint mvMasterMeshID;
    vuint mvSlaveMeshID;

    map<uiint, uiint, less<uiint> >  mmMasterMeshID2Index;
    map<uiint, uiint, less<uiint> >  mmSlaveMeshID2Index;

    vector<CContactNode*> mvConNode;
    map<uiint, uiint, less<uiint> >  mmConNodeID2Index;

    vector<CContactNode*> mvMasterConNode;
    vector<CSkinFace*>    mvFace;         //-------------- マスター面
    map<uiint, uiint, less<uiint> > mmMasterConNodeID2Index;
    map<uiint, uiint, less<uiint> > mmMasterFaceID2Index;

    vector<CContactNode*> mvSlaveConNode;
    vector<CSkinFace*>    mvSlaveFace;    //-------------- スレーブ面
    map<uiint, uiint, less<uiint> > mmSlaveConNodeID2Index;
    map<uiint, uiint, less<uiint> > mmSlaveFaceID2Index;

    // Octree
    CBoundingBox mBoundingBox;
    COctreeKnot moOctreeKnot;
    vector<vector<COctreeKnot*> > mvKnot;


public:
    void setID(const uiint& id) {
        mContactID = id;
    }
    uiint& getID() {
        return mContactID;
    }

    void setLevel(const uiint& level) {
        mLevel= level;
    }
    uiint& getLevel() {
        return mLevel;
    }

    void setProp(const uiint& nProp) {
        mnProp = nProp;
    }
    uiint& getProp() {
        return mnProp;
    }

    void setRank(const uiint& rank) {
        myRank= rank;
    }
    uiint& getRank() {
        return myRank;
    }
    //
    // 伝達率：Film = 線形方程式別の伝達率が入っている
    //
    void setFilm(CFilm* pFilm) {
        mpFilm= pFilm;
    }
    CFilm* getFilm() {
        return mpFilm;
    }
    double& getTransCoeff(const uiint& ieq);


////    void resizeTransmitRank(const uiint& nNumTransmit){ mvTransmitRank.resize(nNumTransmit);}
////    void setTransmitRank(const uiint& i, const uiint& nTransRank){ mvTransmitRank[i]= nTransRank;}
////    uiint& getTransmitRank(const uiint& i){ return mvTransmitRank[i];}
////    vuint& getTransmitRank(){ return mvTransmitRank;}
////    uiint& getNumOfTransmitRank(){ return mNumOfTransmit;}

    void reserveConNode(const uiint& res_size) {
        mvConNode.reserve(res_size);
    }
    void resizeConNode(const uiint& res_size) {
        mvConNode.resize(res_size);
    }
    void addConNode(CContactNode* pConNode, const uiint& id);
    CContactNode* getContactNode(const uiint& index) {
        return mvConNode[index];
    }
    CContactNode* getContactNode_ID(const uiint& id) {
        uiint i=mmConNodeID2Index[id];
        return mvConNode[i];
    }
    uiint getNumOfConNode() {
        return mvConNode.size();
    }

    void addMasterMeshID(const uiint& id);
    uiint& getMasterMeshID(const uiint& index) {
        return mvMasterMeshID[index];
    }
    uiint getNumOfMasterMesh() {
        return mvMasterMeshID.size();
    }
    void addSlaveMeshID(const uiint& id);
    uiint& getSlaveMeshID(const uiint& index) {
        return mvSlaveMeshID[index];
    }
    uiint getNumOfSlaveMesh() {
        return mvSlaveMeshID.size();
    }

    void reserveMasterFace(const uiint& res_size) {
        mvFace.reserve(res_size);
    }
    void resizeMasterFace(const uiint& res_size) {
        mvFace.resize(res_size);
    }
    void addMasterFace(CSkinFace* pFace);
    void addMasterFace(vector<CSkinFace*>& vface);
    void setMasterFace(CSkinFace* pFace, const uiint& index) {
        mvFace[index] = pFace;
    }
    uiint getNumOfMasterFace() {
        return mvFace.size();
    }
    CSkinFace* getMasterFace(const uiint& index) {
        return mvFace[index];
    }
    CSkinFace* getMasterFace_ID(const uiint& id);
    void addMasterConNode(CContactNode* pConNode, const uiint& id);
    void reserveMasterConNode(const uiint& res_size) {
        mvMasterConNode.reserve(res_size);
    }
    void resizeMasterConNode(const uiint& res_size) {
        mvMasterConNode.resize(res_size);
    }
    CContactNode* getMasterConNode(const uiint& index) {
        return mvMasterConNode[index];
    }
    CContactNode* getMasterConNode_ID(const uiint& id) {
        uiint i=mmMasterConNodeID2Index[id];
        return mvMasterConNode[i];
    }
    uiint getNumOfMasterConNode() {
        return mvMasterConNode.size();
    }

    void reserveSlaveFace(const uiint& res_size) {
        mvSlaveFace.reserve(res_size);
    }
    void resizeSlaveFace(const uiint& res_size) {
        mvSlaveFace.resize(res_size);
    }
    void addSlaveFace(CSkinFace* pSlaveFace);
    void addSlaveFace(vector<CSkinFace*>& vface);
    void setSlaveFace(CSkinFace* pSlaveFace, const uiint& index) {
        mvSlaveFace[index]= pSlaveFace;
    }
    uiint getNumOfSlaveFace() {
        return mvSlaveFace.size();
    }
    CSkinFace* getSlaveFace(const uiint& index) {
        return mvSlaveFace[index];
    }
    void addSlaveConNode(CContactNode* pConNode, const uiint& id);
    void reserveSlaveConNode(const uiint& res_size) {
        mvSlaveConNode.reserve(res_size);
    }
    void resizeSlaveConNode(const uiint& res_size) {
        mvSlaveConNode.resize(res_size);
    }
    CContactNode* getSlaveConNode(const uiint& index) {
        return mvSlaveConNode[index];
    }
    CContactNode* getSlaveConNode_ID(const uiint& id) {
        uiint i=mmSlaveConNodeID2Index[id];
        return mvSlaveConNode[i];
    }
    uiint getNumOfSlaveConNode() {
        return mvSlaveConNode.size();
    }


    void setupAggSkinFace();
    void setupEdgeConNode(CContactMesh *pProgConMesh, const uiint& iLevel);
    void setupFaceConNode(CContactMesh *pProgConMesh);
    void setupCoarseConNode(CContactMesh *pProgConMesh);

    void setupSPointOnMFace();
    void setupMPC_Coef();
    uiint getNumOfSlavePoint() {
        return mvSlaveConNode.size();
    }

    void generateOctree(const uiint& maxLayer);

    void deleteProgData();

    void clear();
};
}
#endif
