/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/IndexBucketMesh.h
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
#ifndef _INDEXBUCKET_Mesh_H_a50ffe38_1208_4d65_8605_32ea4cb9cfb1
#define	_INDEXBUCKET_Mesh_H_a50ffe38_1208_4d65_8605_32ea4cb9cfb1
#include "CommonStd.h"
#include "TypeDef.h"
namespace pmw
{
class CIndexBucketMesh
{
public:
    CIndexBucketMesh();
    virtual ~CIndexBucketMesh();
protected:
    uiint maxMeshID, minMeshID;
    uiint defMeshID;
    vuint mvID2MeshIndex;
public:
    virtual void Initialize(const uiint& maxID, const uiint& minID);
    virtual void resizeBucketMesh(const uiint& size);
    virtual void setIndexMesh(const uiint& id, const uiint& index_num);
    virtual uiint& getIndexMesh(const uiint& id) {
        return mvID2MeshIndex[ id - minMeshID];
    }
};
}
#endif	/* _INDEXBUCKET_Mesh_H_ */
