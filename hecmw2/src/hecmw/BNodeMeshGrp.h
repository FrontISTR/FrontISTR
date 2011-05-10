/* 
 * File:   BNodeMeshGrp.h
 * Author: ktakeda
 *
 * BoundaryGroup<BoundaryNodeMesh*> のラッパー
 * # 各階層が全て同一のBNodeMeshGrp* を所有する.
 *
 * Created on 2010/06/25, 14:30
 */
#include "BoundaryGroup.h"
#include "BoundaryNodeMesh.h"

namespace pmw{
#ifndef _BNODEMESH_GRP_H
#define	_BNODEMESH_GRP_H
class CBNodeMeshGrp{
public:
    CBNodeMeshGrp();
    virtual ~CBNodeMeshGrp();

protected:
    BoundaryGroup<CBoundaryNodeMesh*> mGrpBndNodeMesh;

public:
    void reserveBndNodeMesh(const uint& res_size){ mGrpBndNodeMesh.reserve(res_size);}
    void setBndNodeMesh(CBoundaryNodeMesh *pBNodeMesh){ mGrpBndNodeMesh.push(pBNodeMesh);}
    CBoundaryNodeMesh* getBndNodeMeshIX(const uint& index){ return mGrpBndNodeMesh.get_withIndex(index);}
    CBoundaryNodeMesh* getBndNodeMeshID(const uint& id){ return mGrpBndNodeMesh.get_withID(id);}
    uint getNumOfBoundaryNodeMesh(){ return mGrpBndNodeMesh.NumOfBoundary();}
};
#endif	/* _BNODEMESH_GRP_H */
}


