/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   IndexBucketMesh.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _INDEXBUCKET_Mesh_H_a50ffe38_1208_4d65_8605_32ea4cb9cfb1
#define	_INDEXBUCKET_Mesh_H_a50ffe38_1208_4d65_8605_32ea4cb9cfb1
#include "CommonStd.h"
#include "TypeDef.h"
namespace pmw{
class CIndexBucketMesh{
public:
    CIndexBucketMesh();
    virtual ~CIndexBucketMesh();
protected:
    uint maxMeshID, minMeshID;
    uint defMeshID;
    vint mvID2MeshIndex;
public:
    virtual void Initialize(const uint& maxID, const uint& minID);    
    virtual void resizeBucketMesh(const uint& size);                  
    virtual void setIndexMesh(const uint& id, const uint& index_num); 
    virtual int& getIndexMesh(const uint& id){ return mvID2MeshIndex[ id - minMeshID];}
};
}
#endif	/* _INDEXBUCKET_Mesh_H_ */
