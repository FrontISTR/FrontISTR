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
    uint maxNodeID, maxElementID;
    uint minNodeID, minElementID;
    uint defNodeID, defElementID;

    //
    vint mvID2NodeIndex;
    vint mvID2ElementIndex;

public:
    //
    //    virtual void Initialize(const uint& max_nodeID, const uint& min_nodeID,
    //                                const uint& max_elemID, const uint& min_elemID);

    void clearBucketNode();
    void clearBucketElement();

    // resize vector
    void resizeBucketNode(const uint& maxID, const uint& minID);
    void resizeBucketElement(const uint& maxID, const uint& minID);
    void re_resizeBucketNode(const uint& new_maxID);//Nodeが追加された時の再resize

    // setup mvID2Index
    void setIndexNode(const uint& id, const int& index_num);
    void setIndexElement(const uint& id, const int& index_num);

    // get Index
    int& getIndexNode(const uint& id){ return mvID2NodeIndex[id - minNodeID];}
    int& getIndexElement(const uint& id){ return mvID2ElementIndex[id - minElementID];}
};
}
#endif
