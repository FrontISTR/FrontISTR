/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AssyModel.h
|
|                     Written by T.Takeda,    2010/06/01
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
    uint mMaxMeshID,mMinMeshID;
    CIndexBucketMesh moBucketMesh;
    map<uint, uint, less<uint> > mmContactID2Index;
    vector<CMesh*> mvMesh;
    vector<CContactMesh*> mvContactMesh;
    uint mMGLevel;
    CAssyMatrix *mpAssyMatrix;
    CAssyVector *mpAssyVector;
    CAssyVector *mpAssyVector2;
public:
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}
    void setAssyMatrix(CAssyMatrix *assy){ mpAssyMatrix = assy;}
    CAssyMatrix* getAssyMatrix(){ return mpAssyMatrix;}
    void setAssyVector(CAssyVector *vect){ mpAssyVector = vect;}
    CAssyVector* getAssyVector(){ return mpAssyVector;}
    void setAssyVector2(CAssyVector *vect){ mpAssyVector2 = vect;}
    CAssyVector* getAssyVector2(){ return mpAssyVector2;}
    void   resizeMesh(const uint& size){ mvMesh.resize(size);}
    void   setMesh(CMesh *pMesh, const uint& i){  mvMesh[i] = pMesh;}
    CMesh* getMesh(const uint& index){ return mvMesh[index];}
    uint   getNumOfMesh(){ return mvMesh.size();}
    CMesh* getMesh_ID(const uint& id);
    void  addContactMesh(CContactMesh *pContactMesh, const uint& id);
    void  resizeContactMesh(const uint& res_size){  mvContactMesh.resize(res_size);}
    void  setContactMesh(CContactMesh *pContactMesh, const uint& index, const uint& id){ mvContactMesh[index]= pContactMesh; mmContactID2Index[id]=index;}
    CContactMesh* getContactMesh(const uint& index){  return mvContactMesh[index];}
    CContactMesh* getContactMesh_ID(const uint& id){ uint i= mmContactID2Index[id]; return mvContactMesh[i];}
    uint getNumOfContactMesh(){ return mvContactMesh.size();}
    void intializeBucket(const uint& maxID, const uint& minID){ moBucketMesh.Initialize(maxID, minID); }
    void setBucket(const uint& id, const uint& index){ moBucketMesh.setIndexMesh(id, index); }
    void setMaxMeshID(const uint& maxID){ mMaxMeshID= maxID;}
    uint& getMaxMeshID(){ return mMaxMeshID;}
    void setMinMeshID(const uint& minID){ mMinMeshID= minID;}
    uint& getMinMeshID(){ return mMinMeshID;}
};
}
#endif
