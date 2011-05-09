//
//  
//
//
//
//
//			2009.03.24
//			2008.11.10
//			k.Takeda
#ifndef MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6
#define MODEL_HH_AF0F9F8A_E352_4c02_B787_1F8D3BD9D0C6

#include "CommonStd.h"

#include "Mesh.h"
#include "ContactMesh.h"
#include "IndexBucketMesh.h"

namespace pmw{
class CAssyModel
{
public:
    CAssyModel(void);
    virtual ~CAssyModel(void);

protected:
    CIndexBucketMesh moBucketMesh;// Mesh_ID -> index

    vector<CMesh*> mvMesh;
    vector<CContactMesh*> mvContactMesh;

    uint mMGLevel;

public:
    // MultiGrid Level 確認用
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}

    // Mesh
    //void   addMesh(CMesh *pMesh){  mvMesh.push_back(pMesh);}
    void   resizeMesh(const uint& size){ mvMesh.resize(size);}
    void   setMesh(CMesh *pMesh, const uint& i){  mvMesh[i] = pMesh;}
    CMesh* getMesh(const uint& index){ return mvMesh[index];}
    uint   getNumOfMesh(){ return mvMesh.size();}


    // ContactMesh
    void  addContactMesh(CContactMesh *pContactMesh){ mvContactMesh.push_back(pContactMesh);}
    void  resizeContactMesh(const uint& size){  mvContactMesh.resize(size);}
    void  setContactMesh(CContactMesh *pContactMesh, const uint& i){ mvContactMesh[i] = pContactMesh;}

    CContactMesh* getContactMesh(const uint& i){  return mvContactMesh[i];}


    // IndexBucketMesh in AssyModel
    void intializeBucket(const uint& maxID, const uint& minID){ moBucketMesh.Initialize(maxID, minID); }
    void setBucket(const uint& id, const uint& index){ moBucketMesh.setIndexMesh(id, index); }

};
}
#endif
