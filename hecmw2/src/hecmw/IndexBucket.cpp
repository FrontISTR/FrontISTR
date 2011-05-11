//
//  IndexBucket.cpp
//
//		2008.12.15
//		2008.12.12
//		k.Takeda
#include "IndexBucket.h"
using namespace pmw;

#include <iostream>
//
//
CIndexBucket::CIndexBucket()
{
    maxNodeID = 0; maxElementID = 0;
    minNodeID = 0; minElementID = 0;
    defNodeID = 0; defElementID = 0;
}

CIndexBucket::~CIndexBucket()
{
    ////debug
    //std::cout << "~CIndexBucket" << std::endl;
}

// メンバークリア
//
void CIndexBucket::clearBucketNode()
{
    mvID2NodeIndex.clear();
    
    maxNodeID = 0;
    minNodeID = 0;
    defNodeID = 0;
}
void CIndexBucket::clearBucketElement()
{
    mvID2ElementIndex.clear();
    
    maxElementID = 0;
    minElementID = 0;
    defElementID = 0;
}

//// resizeBucket
////
//void CIndexBucket::Initialize(const uint& max_nodeID, const uint& min_nodeID,
//		                      const uint& max_elemID, const uint& min_elemID)
//{
//    maxNodeID = max_nodeID;
//    minNodeID = min_nodeID;
//    maxElementID = max_elemID;
//    minElementID = min_elemID;
//
//    defNodeID = maxNodeID - minNodeID;
//    defElementID = maxElementID - minElementID;
//
//    mvID2NodeIndex.resize(defNodeID + 1);
//    mvID2ElementIndex.resize(defElementID + 1);
//}

// 初期の領域確保
//
void CIndexBucket::resizeBucketNode(const uiint &maxID, const uiint &minID)
{
    maxNodeID = maxID;
    minNodeID = minID;
    defNodeID = maxID - minID;

    mvID2NodeIndex.resize(defNodeID+1);
}
//
// Nodeが追加された時の再resize
//
void CIndexBucket::re_resizeBucketNode(const uiint& new_maxID)
{
    maxNodeID = new_maxID;
    defNodeID = maxNodeID - minNodeID;

    mvID2NodeIndex.resize(defNodeID+1);
}

void CIndexBucket::resizeBucketElement(const uiint &maxID, const uiint &minID)
{
    maxElementID = maxID;
    minElementID = minID;
    defElementID = maxID - minID;

    mvID2ElementIndex.resize(defElementID+1);
}


//
//
void CIndexBucket::setIndexNode(const uiint& id, const int& index_num)
{
    mvID2NodeIndex[id - minNodeID] = index_num;
}

void CIndexBucket::setIndexElement(const uiint &id, const int &index_num)
{
    mvID2ElementIndex[id - minElementID] = index_num;
}


