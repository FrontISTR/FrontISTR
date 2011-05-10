/* 
 * File:   IndexBucketAssyModel.h
 * Author: ktakeda
 *
 * Created on 2009/04/27, 16:49
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
    uint defMeshID;//　<= MaxID - MinID

    vint mvID2MeshIndex;
public:
    virtual void Initialize(const uint& maxID, const uint& minID);    // mvID2Indexの領域確保, 変数初期化
    virtual void resizeBucketMesh(const uint& size);                  // resize vector
    virtual void setIndexMesh(const uint& id, const uint& index_num); // setup mvID2Index
    virtual int& getIndexMesh(const uint& id){ return mvID2MeshIndex[ id - minMeshID];}// get Index
};
}
#endif	/* _INDEXBUCKET_Mesh_H_ */

