//
// IndexBucket.h
//
//
//                          2009.01.06
//                          2008.12.12
//                          k.Takeda

#ifndef INDEX_BUCKET_HH_BE8FA7E_496A_4bc7_9AD0_16E37A0FEEFA
#define INDEX_BUCKET_HH_BE8FA7E_496A_4bc7_9AD0_16E37A0FEEFA

#include "CommonStd.h"
#include "TypeDef.h"

namespace pmw{
class CIndexBucket{
public:
    CIndexBucket();
    virtual ~CIndexBucket();

protected:
    //
    uiint maxNodeID, maxElementID;
    uiint minNodeID, minElementID;
    uiint defNodeID, defElementID;

    //
    vuint mvID2NodeIndex;
    vuint mvID2ElementIndex;

public:
    //
    //    virtual void Initialize(const uint& max_nodeID, const uint& min_nodeID,
    //                                const uint& max_elemID, const uint& min_elemID);

    void clearBucketNode();
    void clearBucketElement();

    // resize vector
    void resizeBucketNode(const uiint& maxID, const uiint& minID);
    void resizeBucketElement(const uiint& maxID, const uiint& minID);
    void re_resizeBucketNode(const uiint& new_maxID);//Nodeが追加された時の再resize

    // setup mvID2Index
    void setIndexNode(const uiint& id, const int& index_num);
    void setIndexElement(const uiint& id, const int& index_num);

    // get Index
    uiint& getIndexNode(const uiint& id){ return mvID2NodeIndex[id - minNodeID];}
    uiint& getIndexElement(const uiint& id){ return mvID2ElementIndex[id - minElementID];}

    // Max Min
    uiint& getMaxNodeID(){ return maxNodeID;}
    uiint& getMinNodeID(){ return minNodeID;}
    uiint& getMaxElementID(){ return maxElementID;}
    uiint& getMinElementID(){ return minElementID;}
};
}
#endif
