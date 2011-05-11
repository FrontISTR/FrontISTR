//
//  AssyModel.h
//
//
//
//
//			2009.10.16
//			2008.11.10
//			k.Takeda
#ifndef MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6
#define MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6

#include "CommonStd.h"
#include "TypeDef.h"

#include "Mesh.h"
#include "ContactMesh.h"
#include "IndexBucketMesh.h"
#include "AssyMatrix.h"

#include <map>

namespace pmw{
class CAssyModel
{
public:
    CAssyModel(void);
    virtual ~CAssyModel(void);

protected:
    uiint mMaxMeshID,mMinMeshID;
    CIndexBucketMesh moBucketMesh;// Mesh_ID -> index
    map<uiint, uiint, less<uiint> > mmContactID2Index;//ContactMesh ID -> index

    vector<CMesh*> mvMesh;
    vector<CContactMesh*> mvContactMesh;

    uiint mMGLevel;

    uiint mNumOfEquation;
    CAssyMatrix **mvAssyMatrix;   // アセンブルマトリックス
    CAssyVector **mvRHSAssyVector;// RHSベクトル
    CAssyVector **mvSolAssyVector;// 解ベクトル

public:
    // MultiGrid Level 確認用
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}

    // Ax=b 管理
    //
    void GeneLinearAlgebra(const vuint& vDOF, CAssyModel *pCoarseAssyModel);//vDOF was each Ax=b DOF
    CAssyMatrix* getAssyMatrix(const uiint& index);   //index is Ax=b number
    CAssyVector* getRHSAssyVector(const uiint& index);//index is Ax=b number
    CAssyVector* getSolutionAssyVector(const uiint& index);//index is Ax=b number
    uiint& getNumOfEquation(){ return mNumOfEquation;}
    
    
    
    // Mesh
    //
    void   resizeMesh(const uiint& size){ mvMesh.resize(size);}
    void   setMesh(CMesh *pMesh, const uiint& i){  mvMesh[i] = pMesh;}
    CMesh* getMesh(const uiint& index){ return mvMesh[index];}
    uiint   getNumOfMesh(){ return mvMesh.size();}
    CMesh* getMesh_ID(const uiint& id);


    // ContactMesh
    //
    void  addContactMesh(CContactMesh *pContactMesh, const uiint& id);
    void  resizeContactMesh(const uiint& res_size){  mvContactMesh.resize(res_size);}
    void  setContactMesh(CContactMesh *pContactMesh, const uiint& index, const uiint& id){ mvContactMesh[index]= pContactMesh; mmContactID2Index[id]=index;}

    CContactMesh* getContactMesh(const uiint& index){  return mvContactMesh[index];}
    CContactMesh* getContactMesh_ID(const uiint& id){ uiint i= mmContactID2Index[id]; return mvContactMesh[i];}
    uiint getNumOfContactMesh(){ return mvContactMesh.size();}

    // IndexBucketMesh in AssyModel
    //
    void intializeBucket(const uiint& maxID, const uiint& minID){ moBucketMesh.Initialize(maxID, minID); }
    void setBucket(const uiint& id, const uiint& index){ moBucketMesh.setIndexMesh(id, index); }

    // MeshIDの最大-最小
    //
    void setMaxMeshID(const uiint& maxID){ mMaxMeshID= maxID;}
    uiint& getMaxMeshID(){ return mMaxMeshID;}
    void setMinMeshID(const uiint& minID){ mMinMeshID= minID;}
    uiint& getMinMeshID(){ return mMinMeshID;}
};
}
#endif
