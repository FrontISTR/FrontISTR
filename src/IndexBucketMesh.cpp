/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   IndexBucketMesh.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
    maxMeshID = 0; minMeshID = 0;
    defMeshID = 0;
}
CIndexBucketMesh::~CIndexBucketMesh()
{
    std::cout << "~CIndexBucketMesh" << std::endl;
}
void CIndexBucketMesh::Initialize(const uint& maxID, const uint& minID)
{
    maxMeshID = maxID;
    minMeshID = minID;
    defMeshID = maxID - minID;    
    mvID2MeshIndex.resize(defMeshID + 1);    
}
void CIndexBucketMesh::resizeBucketMesh(const uint& size)
{
    mvID2MeshIndex.resize(size);
}
void CIndexBucketMesh::setIndexMesh(const uint& id, const uint& index_num)
{
    mvID2MeshIndex[id - minMeshID] = index_num;
}
