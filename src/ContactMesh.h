/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ContactMesh.h
|
|                     Written by T.Takeda,    2010/06/01
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
#include <map>
#include <utility>
typedef std::pair<pmw::CContactNode*, pmw::CContactNode*> PairConNode;
namespace pmw{
class CContactMesh{
public:
    CContactMesh();
    virtual ~CContactMesh();
protected:
    uint mContactID;
    uint mLevel;
    uint myRank;      
    uint transmitRank;
    vuint mvMasterMeshID;
    vuint mvSlaveMeshID; 
    map<uint, uint, less<uint> >  mmMasterMeshID2Index;
    map<uint, uint, less<uint> >  mmSlaveMeshID2Index;
    vector<CContactNode*> mvConNode;
    map<uint, uint, less<uint> >  mmConNodeID2Index;
    vector<CContactNode*> mvMasterConNode;
    vector<CSkinFace*>    mvFace;         
    map<uint, uint, less<uint> > mmMasterConNodeID2Index;
    map<uint, uint, less<uint> > mmMasterFaceID2Index;
    vector<CContactNode*> mvSlaveConNode;
    vector<CSkinFace*>    mvSlaveFace;   
    map<uint, uint, less<uint> > mmSlaveConNodeID2Index;
    map<uint, uint, less<uint> > mmSlaveFaceID2Index;
    CBoundingBox mBoundingBox;
    COctreeKnot moOctreeKnot;            
    vector<vector<COctreeKnot*> > mvKnot;
public:
    void setID(const uint& id){ mContactID = id;}
    uint& getID(){ return mContactID;}
    void setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
    void setRank(const uint& rank){ myRank= rank;}
    uint& getRank(){ return myRank;}
    void setTransmitRank(const uint& rank){ transmitRank= rank;}
    uint& getTransmitRank(){ return transmitRank;}
    void reserveConNode(const uint& res_size){ mvConNode.reserve(res_size);}
    void resizeConNode(const uint& res_size){ mvConNode.resize(res_size);}
    void addConNode(CContactNode* pConNode, const uint& id);
    CContactNode* getContactNode(const uint& index){ return mvConNode[index];}
    CContactNode* getContactNode_ID(const uint& id){ uint i=mmConNodeID2Index[id]; return mvConNode[i];}
    uint getNumOfConNode(){ return mvConNode.size();}
    void addMasterMeshID(const uint& id);
    uint& getMasterMeshID(const uint& index){return mvMasterMeshID[index];}
    uint getNumOfMasterMesh(){ return mvMasterMeshID.size();}
    void addSlaveMeshID(const uint& id);
    uint& getSlaveMeshID(const uint& index){ return mvSlaveMeshID[index];}
    uint getNumOfSlaveMesh(){ return mvSlaveMeshID.size();}
    void reserveMasterFace(const uint& res_size){ mvFace.reserve(res_size);}
    void resizeMasterFace(const uint& res_size){ mvFace.resize(res_size);}
    void addMasterFace(CSkinFace* pFace);
    void addMasterFace(vector<CSkinFace*>& vface);
    void setMasterFace(CSkinFace* pFace, const uint& index){ mvFace[index] = pFace;}
    uint getNumOfMasterFace(){ return mvFace.size();}
    CSkinFace* getMasterFace(const uint& index){ return mvFace[index];}
    CSkinFace* getMasterFace_ID(const uint& id);
    void addMasterConNode(CContactNode* pConNode, const uint& id);
    void reserveMasterConNode(const uint& res_size){ mvMasterConNode.reserve(res_size);}
    void resizeMasterConNode(const uint& res_size){ mvMasterConNode.resize(res_size);}
    CContactNode* getMasterConNode(const uint& index){ return mvMasterConNode[index];}
    CContactNode* getMasterConNode_ID(const uint& id){ uint i=mmMasterConNodeID2Index[id]; return mvMasterConNode[i];}
    uint getNumOfMasterConNode(){ return mvMasterConNode.size();}
    void reserveSlaveFace(const uint& res_size){ mvSlaveFace.reserve(res_size);}
    void resizeSlaveFace(const uint& res_size){ mvSlaveFace.resize(res_size);}
    void addSlaveFace(CSkinFace* pSlaveFace);
    void addSlaveFace(vector<CSkinFace*>& vface);
    void setSlaveFace(CSkinFace* pSlaveFace, const uint& index){ mvSlaveFace[index]= pSlaveFace;}
    uint getNumOfSlaveFace(){ return mvSlaveFace.size();}
    CSkinFace* getSlaveFace(const uint& index){ return mvSlaveFace[index];}
    void addSlaveConNode(CContactNode* pConNode, const uint& id);
    void reserveSlaveConNode(const uint& res_size){ mvSlaveConNode.reserve(res_size);}
    void resizeSlaveConNode(const uint& res_size){ mvSlaveConNode.resize(res_size);}
    CContactNode* getSlaveConNode(const uint& index){ return mvSlaveConNode[index];}
    CContactNode* getSlaveConNode_ID(const uint& id){ uint i=mmSlaveConNodeID2Index[id]; return mvSlaveConNode[i];}
    uint getNumOfSlaveConNode(){ return mvSlaveConNode.size();}
    void setupAggSkinFace();
    void setupEdgeConNode(CContactMesh *pProgConMesh);
    void setupFaceConNode(CContactMesh *pProgConMesh);
    void setupCoarseConNode(CContactMesh *pProgConMesh);
    void setupSPointOnMFace();
    void setupMPC_Coef();     
    uint getNumOfSlavePoint(){ return mvSlaveConNode.size();}
    void generateOctree(const uint& maxLayer);
};
}
#endif
