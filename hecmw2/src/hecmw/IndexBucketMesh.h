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
    uiint maxMeshID, minMeshID;
    uiint defMeshID;//　<= MaxID - MinID

    vuint mvID2MeshIndex;
public:
    virtual void Initialize(const uiint& maxID, const uiint& minID);    // mvID2Indexの領域確保, 変数初期化
    virtual void resizeBucketMesh(const uiint& size);                  // resize vector
    virtual void setIndexMesh(const uiint& id, const uiint& index_num); // setup mvID2Index
    virtual uiint& getIndexMesh(const uiint& id){ return mvID2MeshIndex[ id - minMeshID];}// get Index
};
}
#endif	/* _INDEXBUCKET_Mesh_H_ */

