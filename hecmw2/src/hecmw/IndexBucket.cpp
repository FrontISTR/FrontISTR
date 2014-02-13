/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/IndexBucket.cpp
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
#include "IndexBucket.h"
using namespace pmw;
#include <iostream>
CIndexBucket::CIndexBucket()
{
    maxNodeID = 0;
    maxElementID = 0;
    minNodeID = 0;
    minElementID = 0;
    defNodeID = 0;
    defElementID = 0;
}
CIndexBucket::~CIndexBucket()
{
}
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
void CIndexBucket::resizeBucketNode(const uiint &maxID, const uiint &minID)
{
    maxNodeID = maxID;
    minNodeID = minID;
    defNodeID = maxID - minID;
    mvID2NodeIndex.resize(defNodeID+1);
}
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
void CIndexBucket::setIndexNode(const uiint& id, const int& index_num)
{
    mvID2NodeIndex[id - minNodeID] = index_num;
}
void CIndexBucket::setIndexElement(const uiint &id, const int &index_num)
{
    mvID2ElementIndex[id - minElementID] = index_num;
}
