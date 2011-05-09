/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   IndexBucket.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "IndexBucket.h"
using namespace pmw;
#include <iostream>
CIndexBucket::CIndexBucket()
{
    maxNodeID = 0; maxElementID = 0;
    minNodeID = 0; minElementID = 0;
    defNodeID = 0; defElementID = 0;
}
CIndexBucket::~CIndexBucket()
{
    std::cout << "~CIndexBucket" << std::endl;
}
void CIndexBucket::Initialize(const uint& max_nodeID, const uint& min_nodeID,
		                      const uint& max_elemID, const uint& min_elemID)
{
    maxNodeID = max_nodeID;
    minNodeID = min_nodeID;
    maxElementID = max_elemID;
    minElementID = min_elemID;
    defNodeID = maxNodeID - minNodeID;
    defElementID = maxElementID - minElementID;
    mvID2NodeIndex.resize(defNodeID + 1);
    mvID2ElementIndex.resize(defElementID + 1);
}
void CIndexBucket::resizeBucketNode(const uint &maxID, const uint &minID)
{
    maxNodeID = maxID;
    minNodeID = minID;
    defNodeID = maxID - minID;
    mvID2NodeIndex.resize(defNodeID+1);
}
void CIndexBucket::resizeBucketElement(const uint &maxID, const uint &minID)
{
    maxElementID = maxID;
    minElementID = minID;
    defElementID = maxID - minID;
    mvID2ElementIndex.resize(defElementID+1);
}
void CIndexBucket::setIndexNode(const uint& id, const int& index_num)
{
    mvID2NodeIndex[id - minNodeID] = index_num;
}
void CIndexBucket::setIndexElement(const uint &id, const int &index_num)
{
    mvID2ElementIndex[id - minElementID] = index_num;
}
