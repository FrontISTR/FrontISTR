/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/AssyModel.h
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
#ifndef MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6
#define MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6
#include "CommonStd.h"
#include "TypeDef.h"
#include "Mesh.h"
#include "ContactMesh.h"
#include "IndexBucketMesh.h"
#include "AssyMatrix.h"
#include <map>
#include <utility>
namespace pmw
{
class CAssyModel
{
public:
    CAssyModel(void);
    virtual ~CAssyModel(void);
protected:
    uiint mMaxMeshID,mMinMeshID;

    CIndexBucketMesh moBucketMesh;

    map<uiint, uiint, less<uiint> > mmContactID2Index;

    vector<CMesh*> mvMesh;
    vector<CContactMesh*> mvContactMesh;
    uiint mMGLevel;

    uiint mnNumOfGlobalComm;              // 全領域のCommMesh2数
    vector<pair<uiint,uiint> > mvRankPair;// 全領域のCommMesh2のmyRank,transRank : [テーブルNum].first=rank0, [テーブルNum].second=rank1
    vuint mvMeshID_CommID;                // グローバルCommIDごとのMeshID番号


    uiint mNumOfEquation;// Number of Linear_Algebra_Equation

    CAssyMatrix **mvAssyMatrix;
    CAssyVector **mvRHSAssyVector;
    CAssyVector **mvSolAssyVector;

public:
    void setMGLevel(const uiint& mgLevel) {
        mMGLevel= mgLevel;
    }
    uiint& getMGLevel() {
        return mMGLevel;
    }

    void GeneLinearAlgebra(vvuint& vvDOF, CAssyModel *pCoarseAssyModel);//---- 方程式Ax=b生成

    CAssyMatrix* getAssyMatrix(const uiint& index);
    CAssyVector* getRHSAssyVector(const uiint& index);
    CAssyVector* getSolutionAssyVector(const uiint& index);


    uiint& getNumOfEquation() {
        return mNumOfEquation;
    }

    void   resizeMesh(const uiint& size) {
        mvMesh.resize(size);
    }
    void   setMesh(CMesh *pMesh, const uiint& i) {
        mvMesh[i] = pMesh;
    }
    CMesh* getMesh(const uiint& index) {
        return mvMesh[index];
    }
    uiint   getNumOfMesh() {
        return mvMesh.size();
    }
    CMesh* getMesh_ID(const uiint& id);

    uiint getIndex_of_Mesh(const uiint& id) {
        return moBucketMesh.getIndexMesh(id);   //--Mesh ID番号から、Mesh index番号
    }
    uiint getID_of_Mesh(const uiint& index) {
        return mvMesh[index]->getMeshID();   //--Mesh index番号から、Mesh ID番号
    }

    //--
    // 接合面:ContactMesh
    //--
    void  addContactMesh(CContactMesh *pContactMesh, const uiint& id);
    void  resizeContactMesh(const uiint& res_size) {
        mvContactMesh.resize(res_size);
    }
    void  setContactMesh(CContactMesh *pContactMesh, const uiint& index, const uiint& id) {
        mvContactMesh[index]= pContactMesh;
        mmContactID2Index[id]=index;
    }
    CContactMesh* getContactMesh(const uiint& index) {
        return mvContactMesh[index];
    }
    CContactMesh* getContactMesh_ID(const uiint& id) {
        uiint i= mmContactID2Index[id];
        return mvContactMesh[i];
    }
    uiint getNumOfContactMesh() {
        return mvContactMesh.size();
    }
    uiint getContactID(const uiint& index) {
        return mvContactMesh[index]->getID();
    }

    void intializeBucket(const uiint& maxID, const uiint& minID) {
        moBucketMesh.Initialize(maxID, minID);
    }
    void setBucket(const uiint& id, const uiint& index) {
        moBucketMesh.setIndexMesh(id, index);
    }


    void setMaxMeshID(const uiint& maxID) {
        mMaxMeshID= maxID;
    }
    uiint& getMaxMeshID() {
        return mMaxMeshID;
    }

    void setMinMeshID(const uiint& minID) {
        mMinMeshID= minID;
    }
    uiint& getMinMeshID() {
        return mMinMeshID;
    }

    bool isSelfMesh(const uiint& id);// IDのMeshは存在するか？

    //Global通信テーブル数, 通信テーブルのrankペア(互いに通信する相手)
    uiint& getNumOfGlobalCommMesh2() {
        return mnNumOfGlobalComm;
    }
    uiint& getGlobalPairRank_1st(const uiint& nCommID) {
        return mvRankPair[nCommID].first;
    }
    uiint& getGlobalPairRank_2nd(const uiint& nCommID) {
        return mvRankPair[nCommID].second;
    }
    //Global通信テーブルは、どのパーツ(Mesh)に付属しているか.
    uiint& getMeshID_with_CommID(const uiint& nCommID) {
        return mvMeshID_CommID[nCommID];
    }

    void setNumGlobalCommMesh2(const uiint& nNumOfGlobalComm);
    void setGlobalPairRank(const uiint& nCommID, const uiint& nFirstRank, const uiint& nSecondRank);
    void setMeshID_with_CommID(const uiint& nCommID, const uiint& nMeshID);
};
}
#endif
