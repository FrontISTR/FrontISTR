/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/IndexBucket.h
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
#ifndef INDEX_BUCKET_HH_BE8FA7E_496A_4bc7_9AD0_16E37A0FEEFA
#define INDEX_BUCKET_HH_BE8FA7E_496A_4bc7_9AD0_16E37A0FEEFA
#include "CommonStd.h"
#include "TypeDef.h"
namespace pmw
{
class CIndexBucket
{
public:
    CIndexBucket();
    virtual ~CIndexBucket();
protected:
    uiint maxNodeID, maxElementID;
    uiint minNodeID, minElementID;
    uiint defNodeID, defElementID;
    vuint mvID2NodeIndex;
    vuint mvID2ElementIndex;
public:
    void clearBucketNode();
    void clearBucketElement();
    void resizeBucketNode(const uiint& maxID, const uiint& minID);
    void resizeBucketElement(const uiint& maxID, const uiint& minID);
    void re_resizeBucketNode(const uiint& new_maxID);
    void setIndexNode(const uiint& id, const int& index_num);
    void setIndexElement(const uiint& id, const int& index_num);
    uiint& getIndexNode(const uiint& id) {
        return mvID2NodeIndex[id - minNodeID];
    }
    uiint& getIndexElement(const uiint& id) {
        return mvID2ElementIndex[id - minElementID];
    }
    uiint& getMaxNodeID() {
        return maxNodeID;
    }
    uiint& getMinNodeID() {
        return minNodeID;
    }
    uiint& getMaxElementID() {
        return maxElementID;
    }
    uiint& getMinElementID() {
        return minElementID;
    }
};
}
#endif
