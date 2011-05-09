//
//
//              2009.4.27
//              2009.4.27
//              k.Takeda

#include "IndexBucketMesh.h"
using namespace pmw;

#include <iostream>
// コンストラクター&デストラクター
//
CIndexBucketMesh::CIndexBucketMesh()
{
    maxMeshID = 0; minMeshID = 0;
    defMeshID = 0;
}

CIndexBucketMesh::~CIndexBucketMesh()
{
    //debug
    std::cout << "~CIndexBucketMesh" << std::endl;
}

// 初期化::mvID2Indexの領域確保,変数初期化
//
void CIndexBucketMesh::Initialize(const uint& maxID, const uint& minID)
{
    maxMeshID = maxID;
    minMeshID = minID;
    
    defMeshID = maxID - minID;    
    mvID2MeshIndex.resize(defMeshID + 1);    
}

// mvID2Index の resize
//
void CIndexBucketMesh::resizeBucketMesh(const uint& size)
{
    mvID2MeshIndex.resize(size);
}

// ID -> Index のセットアップ
//
void CIndexBucketMesh::setIndexMesh(const uint& id, const uint& index_num)
{
    mvID2MeshIndex[id - minMeshID] = index_num;
}


