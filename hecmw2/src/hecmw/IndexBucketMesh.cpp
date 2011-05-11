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
    ////debug
    //std::cout << "~CIndexBucketMesh" << std::endl;
}

// 初期化::mvID2Indexの領域確保,変数初期化
//
void CIndexBucketMesh::Initialize(const uiint& maxID, const uiint& minID)
{
    maxMeshID = maxID;
    minMeshID = minID;
    
    defMeshID = maxID - minID;    
    mvID2MeshIndex.resize(defMeshID + 1);    
}

// mvID2Index の resize
//
void CIndexBucketMesh::resizeBucketMesh(const uiint& size)
{
    mvID2MeshIndex.resize(size);
}

// ID -> Index のセットアップ
//
void CIndexBucketMesh::setIndexMesh(const uiint& id, const uiint& index_num)
{
    mvID2MeshIndex[id - minMeshID] = index_num;
}


