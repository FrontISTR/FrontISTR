/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AssyModel.h
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
    CIndexBucketMesh moBucketMesh;
    map<uiint, uiint, less<uiint> > mmContactID2Index;
    vector<CMesh*> mvMesh;
    vector<CContactMesh*> mvContactMesh;
    uiint mMGLevel;
    uiint mNumOfEquation;
    CAssyMatrix **mvAssyMatrix;   
    CAssyVector **mvRHSAssyVector;
    CAssyVector **mvSolAssyVector;
public:
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}
    void GeneLinearAlgebra(const vuint& vDOF, CAssyModel *pCoarseAssyModel);
    CAssyMatrix* getAssyMatrix(const uiint& index);   
    CAssyVector* getRHSAssyVector(const uiint& index);
    CAssyVector* getSolutionAssyVector(const uiint& index);
    uiint& getNumOfEquation(){ return mNumOfEquation;}
    void   resizeMesh(const uiint& size){ mvMesh.resize(size);}
    void   setMesh(CMesh *pMesh, const uiint& i){  mvMesh[i] = pMesh;}
    CMesh* getMesh(const uiint& index){ return mvMesh[index];}
    uiint   getNumOfMesh(){ return mvMesh.size();}
    CMesh* getMesh_ID(const uiint& id);
    void  addContactMesh(CContactMesh *pContactMesh, const uiint& id);
    void  resizeContactMesh(const uiint& res_size){  mvContactMesh.resize(res_size);}
    void  setContactMesh(CContactMesh *pContactMesh, const uiint& index, const uiint& id){ mvContactMesh[index]= pContactMesh; mmContactID2Index[id]=index;}
    CContactMesh* getContactMesh(const uiint& index){  return mvContactMesh[index];}
    CContactMesh* getContactMesh_ID(const uiint& id){ uiint i= mmContactID2Index[id]; return mvContactMesh[i];}
    uiint getNumOfContactMesh(){ return mvContactMesh.size();}
    void intializeBucket(const uiint& maxID, const uiint& minID){ moBucketMesh.Initialize(maxID, minID); }
    void setBucket(const uiint& id, const uiint& index){ moBucketMesh.setIndexMesh(id, index); }
    void setMaxMeshID(const uiint& maxID){ mMaxMeshID= maxID;}
    uiint& getMaxMeshID(){ return mMaxMeshID;}
    void setMinMeshID(const uiint& minID){ mMinMeshID= minID;}
    uiint& getMinMeshID(){ return mMinMeshID;}
};
}
#endif
