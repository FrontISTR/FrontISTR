/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/IndexBucketMesh.cpp
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
#include "IndexBucketMesh.h"
using namespace pmw;
#include <iostream>
CIndexBucketMesh::CIndexBucketMesh()
{
    maxMeshID = 0;
    minMeshID = 0;
    defMeshID = 0;
}
CIndexBucketMesh::~CIndexBucketMesh()
{
}
void CIndexBucketMesh::Initialize(const uiint& maxID, const uiint& minID)
{
    maxMeshID = maxID;
    minMeshID = minID;
    defMeshID = maxID - minID;
    mvID2MeshIndex.resize(defMeshID + 1);
}
void CIndexBucketMesh::resizeBucketMesh(const uiint& size)
{
    mvID2MeshIndex.resize(size);
}
void CIndexBucketMesh::setIndexMesh(const uiint& id, const uiint& index_num)
{
    mvID2MeshIndex[id - minMeshID] = index_num;
}
